import os
import shutil
import time
from FragmentKnitwork.Quilter.align import (alignment, placement,
                                            reverse_alignment)
from FragmentKnitwork.utils import quilterConfig as config
from FragmentKnitwork.utils.quilterUtils import (create_alignment_dirs,
                                                 get_ref_substructures,
                                                 move_files, add_props_from_dict)
from FragmentKnitwork.utils.fragmentUtils import filt_by_bool
from FragmentKnitwork.utils.utils import dump_json, load_json, order_by_lst
from FragmentKnitwork.Quilter.checkPlacedMol import process_molecule
from joblib import Parallel, delayed
from rdkit import Chem


class Pipeline:

    def __init__(self,
                 smiles,
                 fragmentA,
                 fragmentB,
                 target,
                 pdb_file,
                 cmpd_ids,
                 ref_subnodes,
                 ref_synthons,
                 output_dir,
                 working_dir,
                 scoring_dir,
                 merge_type,
                 substructure_dir,
                 similarities=None,
                 used_subnodes=None,
                 used_synthons=None,
                 original_similarities=None,
                 original_names=None,
                 rev_synthons=None,
                 intermed_nodes1=None,
                 intermed_nodes2=None,
                 limit_num_run=None,
                 n_cpus=config.N_CPUS,
                 move_files=None,
                 ):
        # for all merge types
        self.smiles = smiles
        self.mols = [Chem.MolFromSmiles(smi) for smi in self.smiles]
        self.fragmentA = fragmentA
        self.fragmentB = fragmentB
        self.target = target
        self.pdb_file = pdb_file
        self.cmpd_ids = cmpd_ids
        self.ref_subnodes = ref_subnodes
        self.ref_synthons = ref_synthons

        # for impure merges
        self.similarities = similarities
        self.used_subnodes = used_subnodes
        self.used_synthons = used_synthons

        # for reverse merges
        self.original_similarities = original_similarities
        self.original_names = original_names
        self.rev_synthons = rev_synthons
        self.intermed_nodes1 = intermed_nodes1
        self.intermed_nodes2 = intermed_nodes2

        # for all merge types
        self.output_dir = output_dir
        self.working_dir = working_dir
        self.scoring_dir = scoring_dir
        self.merge_type = merge_type
        self.substructure_dir = substructure_dir

        # threshold options
        self.limit_num_run = limit_num_run
        self.n_cpus = n_cpus
        self.move_files = move_files

        # create unique names for merges
        self.names = None
        self.get_names()

        # load the reference substructures to use for placement
        self.ref_subnode_mols, self.ref_synthon_mols = [], []
        self.get_ref_subnode_mols()

        # limit the number of mols for filtering (if requested)
        if self.limit_num_run:
            print('Limiting number of molecules for filtering')
            self.filter_by_limit()
            print(len(self.mols), 'mols for filtering')

        # create dirs for pair results
        self.pair, self.pair_working_dir, self.pair_output_dir = create_alignment_dirs(self.fragmentA, self.fragmentB,
                                                                                       self.working_dir, self.output_dir)
        if os.path.exists(os.path.join(self.pair_output_dir, 'passing_data.json')):
            print(f'Pair {self.pair} already run')

        # dump input data
        self.dump_input_data()

    def get_names(self):
        N = len(self.smiles)
        # if self.merge_type != 'reverse_merge':
        self.names = [f"{self.fragmentA}-{self.fragmentB}_{idx}" for idx in range(N)]

    def filter_by_limit(self):
        """
        If requested, limit number of mols to run through filtering. Done by ordered similarities of replace substructure
        if a bioisosteric merge.

        :return:
        """
        if self.merge_type == 'pure_merge':
            self.names = self.names[:self.limit_num_run]
            self.smiles = self.smiles[:self.limit_num_run]
            self.mols = self.mols[:self.limit_num_run]
            self.cmpd_ids = self.cmpd_ids[:self.limit_num_run]
            self.ref_subnodes = self.ref_subnodes[:self.limit_num_run]
            self.ref_synthons = self.ref_synthons[:self.limit_num_run]
            self.ref_subnode_mols = self.ref_subnode_mols[:self.limit_num_run]
            self.ref_synthon_mols = self.ref_synthon_mols[:self.limit_num_run]

        if self.merge_type == 'impure_merge' or self.merge_type == 'reverse_merge':
            self.names = self.names[:self.limit_num_run]  # doesn't need re-ordering
            self.smiles = order_by_lst(self.smiles, self.similarities)[:self.limit_num_run]
            self.mols = order_by_lst(self.mols, self.similarities)[:self.limit_num_run]
            self.cmpd_ids = order_by_lst(self.cmpd_ids, self.similarities)[:self.limit_num_run]
            self.ref_subnodes = order_by_lst(self.ref_subnodes, self.similarities)[:self.limit_num_run]
            self.ref_synthons = order_by_lst(self.ref_synthons, self.similarities)[:self.limit_num_run]
            self.ref_subnode_mols = order_by_lst(self.ref_subnode_mols, self.similarities)[:self.limit_num_run]
            self.ref_synthon_mols = order_by_lst(self.ref_synthon_mols, self.similarities)[:self.limit_num_run]
            self.used_subnodes = order_by_lst(self.used_subnodes, self.similarities)[:self.limit_num_run]
            self.used_synthons = order_by_lst(self.used_synthons, self.similarities)[:self.limit_num_run]

            if self.merge_type == 'reverse_merge':
                self.original_names = order_by_lst(self.original_names, self.similarities)[:self.limit_num_run]
                self.original_smiles = order_by_lst(self.original_smiles, self.similarities)[:self.limit_num_run]
                self.original_similarities = order_by_lst(self.original_similarities, self.similarities)[:self.limit_num_run]
                self.rev_synthons = order_by_lst(self.rev_synthons, self.similarities)[:self.limit_num_run]
                self.intermed_nodes1 = order_by_lst(self.intermed_nodes1, self.similarities)[:self.limit_num_run]
                if self.intermed_nodes2:
                    self.intermed_nodes2 = order_by_lst(self.intermed_nodes2, self.similarities)[:self.limit_num_run]

            self.similarities.sort(reverse=True)
            self.similarities = self.similarities[:self.limit_num_run]

    def get_data_dict(self):
        """

        :return:
        """
        data = {'names': self.names,
                'smiles': self.smiles,
                'cmpd_ids': self.cmpd_ids,
                'ref_subnodes': self.ref_subnodes,
                'ref_synthons': self.ref_synthons}

        if self.merge_type == 'impure_merge' or self.merge_type == 'reverse_merge':
            data.update({'similarities': self.similarities,
                        'used_subnodes': self.used_subnodes,
                        'used_synthons': self.used_synthons})

        if self.merge_type == 'reverse_merge':
            data.update({'original_smiles': self.original_smiles,
                        'original_names': self.original_names,
                        'original_similarities': self.original_similarities,
                        'rev_synthons': self.rev_synthons,
                        'intermed_nodes1': self.intermed_nodes1})
            if self.intermed_nodes2:
                data.update({'intermed_nodes2': self.intermed_nodes2})

        return data

    def dump_input_data(self):
        """

        :return:
        """
        self.data = self.get_data_dict()

        # write the input data
        dump_json(self.data, os.path.join(self.pair_output_dir, 'input_data.json'))
        print('Input data dumped')

    def get_ref_subnode_mols(self):
        """
        Load the reference substructures to use for the alignment (there may be multiple possible matches to try). The
        ref substructures should be saved as files in the dir created by runEnumeration.py

        :return:
        """
        for subnode, synthon in zip(self.ref_subnodes, self.ref_synthons):

            if self.merge_type != 'reverse_merge':
                ref_subnode = get_ref_substructures(subnode, self.fragmentA, substructure_dir=self.substructure_dir, isSynthon=False)
                self.ref_subnode_mols.append(ref_subnode)
            if self.merge_type == 'reverse_merge':
                ref_subnode = get_ref_substructures(subnode, self.fragmentA, substructure_dir=self.substructure_dir, isSynthon=True)
                self.ref_subnode_mols.extend(ref_subnode)

            ref_synthon = get_ref_substructures(synthon, self.fragmentB, substructure_dir=self.substructure_dir, isSynthon=True)
            self.ref_synthon_mols.append(ref_synthon)

    def run_filtering(self):
        """
        Run the relevant filtering/conformer generation method depending on the type of merge

        :return:
        """
        start = time.time()
        if self.merge_type == 'pure_merge':
            self.filter_results = Parallel(n_jobs=self.n_cpus, backend='multiprocessing')(
                delayed(placement)(
                    smi, name, ref_subnode, ref_synthon, self.pdb_file, self.pair_working_dir
                ) for smi, name, ref_subnode, ref_synthon, in
                zip(self.smiles, self.names, self.ref_subnode_mols, self.ref_synthon_mols)
            )

        if self.merge_type == 'impure_merge':
            self.filter_results = Parallel(n_jobs=self.n_cpus, backend='multiprocessing')(
                delayed(alignment)(
                    smi, name, merge, ref_subnode, ref_synthon, used_subnode, used_synthon,
                    self.pdb_file, self.pair_working_dir
                ) for smi, name, merge, ref_subnode, ref_synthon, used_subnode, used_synthon in
                zip(self.smiles, self.names, self.mols, self.ref_subnode_mols, self.ref_synthon_mols, self.used_subnodes, self.used_synthons)
            )

        if self.merge_type == 'reverse_merge':
            if not self.intermed_nodes2:
                self.filter_results = Parallel(n_jobs=self.n_cpus, backend='multiprocessing')(
                    delayed(reverse_alignment)(
                        smi, original_smi, intermed_smi1, name, merge, ref_subnode, ref_synthon, rev_synthon,
                        used_synthon, self.pdb_file, self.pair_working_dir
                    ) for smi, original_smi, intermed_smi1, name, merge, ref_subnode, ref_synthon, rev_synthon, used_synthon in
                    zip(self.smiles, self.original_smiles, self.intermed_nodes1, self.names, self.mols, self.ref_subnode_mols, self.ref_synthon_mols,
                        self.rev_synthons, self.used_synthons)
                )
            if self.intermed_nodes2:
                self.filter_results = Parallel(n_jobs=self.n_cpus, backend='multiprocessing')(
                    delayed(reverse_alignment)(
                        smi, original_smi, intermed_smi1, name, merge, ref_subnode, ref_synthon, rev_synthon,
                        used_synthon, self.pdb_file, self.pair_working_dir, intermed_smi2
                    ) for smi, original_smi, intermed_smi1, name, merge, ref_subnode, ref_synthon, rev_synthon, used_synthon, intermed_smi2 in
                    zip(self.smiles, self.original_smiles, self.intermed_nodes1, self.names, self.mols, self.ref_subnode_mols, self.ref_synthon_mols,
                        self.rev_synthons, self.used_synthons, self.intermed_nodes2)
                )
        end = time.time()
        self.total_filter_time = round(end - start, 2)
        self.process_filtered_results()

    def save_filtered_results(self):
        """

        :return:
        """
        # results contain: passes, aligned_mols, results_dirs, errors, timings
        print('Processing filter results')
        self.filter_passes, self.filter_mols, self.results_dirs, self.errors, self.timings = [], [], [], [], []
        for res in self.filter_results:
            self.filter_passes.append(res[0])
            self.filter_mols.append(res[1])
            self.results_dirs.append(res[2])
            self.errors.append(res[3])
            self.timings.append(res[4])

        # write summary of results to file
        self.data.update({'passes': self.filter_passes,
                          'errors': self.errors,
                          'timings': self.timings,
                          'total_filter_time': self.total_filter_time,
                          'n_cpus': self.n_cpus})

        dump_json(self.data, os.path.join(self.pair_output_dir, 'passing_data.json'))

        print('Moving successful results files to output dir')
        self.move_results_files()

    def process_filtered_results(self):
        """

        :return:
        """
        # save the results to output json file (and basic processing)
        self.save_filtered_results()

        # filter according to the filter results
        self.names = filt_by_bool(self.names, self.filter_passes)
        self.mols = filt_by_bool(self.mols, self.filter_passes)
        self.smiles = filt_by_bool(self.smiles, self.filter_passes)
        self.cmpd_ids = filt_by_bool(self.cmpd_ids, self.filter_passes)
        self.ref_subnodes = filt_by_bool(self.ref_subnodes, self.filter_passes)
        self.ref_synthons = filt_by_bool(self.ref_synthons, self.filter_passes)
        self.ref_subnode_mols = filt_by_bool(self.ref_subnode_mols, self.filter_passes)
        self.ref_synthon_mols = filt_by_bool(self.ref_synthon_mols, self.filter_passes)
        self.results_dirs = filt_by_bool(self.results_dirs, self.filter_passes)
        self.filter_mols = filt_by_bool(self.filter_mols, self.filter_passes)

        if self.merge_type == 'impure_merge' or self.merge_type == 'reverse_merge':
            self.used_subnodes = filt_by_bool(self.used_subnodes, self.filter_passes)
            self.used_synthons = filt_by_bool(self.used_synthons, self.filter_passes)

            if self.merge_type == 'reverse_merge':
                self.original_names = filt_by_bool(self.original_names, self.filter_passes)
                self.original_smiles = filt_by_bool(self.original_smiles, self.filter_passes)
                self.original_similarities = filt_by_bool(self.original_similarities, self.filter_passes)
                self.rev_synthons = filt_by_bool(self.rev_synthons, self.filter_passes)
                self.intermed_nodes1 = filt_by_bool(self.intermed_nodes1, self.filter_passes)

                if self.intermed_nodes2:
                    self.intermed_nodes2 = filt_by_bool(self.intermed_nodes2, self.filter_passes)

            self.similarities = filt_by_bool(self.similarities, self.filter_passes)

    def move_results_files(self):
        """

        :return:
        """
        if self.move_files:
            results_dirs = [i for i in self.results_dirs if i]
            move_files(results_dirs, self.pair_output_dir, True)
            shutil.rmtree(self.pair_working_dir)

    def get_files_for_scoring(self):
        """

        :return:
        """
        self.holo_files, self.json_files, self.interaction_files, self.mol_files = [], [], [], []

        for results_dir, name in zip(self.results_dirs, self.names):
            name_hyph = name.replace('_', '-')
            holo_file = os.path.join(self.pair_output_dir, name_hyph, f"{name_hyph}.holo_minimised.pdb")
            json_file = os.path.join(self.pair_output_dir, name_hyph, f"{name_hyph}.minimised.json")
            interaction_file = os.path.join(self.pair_output_dir, name_hyph, f"{name_hyph}_interactions.json")
            mol_file = os.path.join(self.pair_output_dir, name_hyph, f"{name_hyph}.minimised.mol")
            self.holo_files.append(holo_file)
            self.json_files.append(json_file)
            self.interaction_files.append(interaction_file)
            self.mol_files.append(mol_file)

        self.fragment_ints_data = load_json(os.path.join(self.substructure_dir, 'fragment_interactions.json'))

    def run_scoring(self):
        """

        :return:
        """
        # get dict to store scoring data
        self.scoring_data = self.get_data_dict()

        # get the files needed for scoring
        self.get_files_for_scoring()

        # run scoring
        print('Running scoring of molecules')
        self.scoring_results = Parallel(n_jobs=self.n_cpus, backend='multiprocessing')(
            delayed(process_molecule)(
                mol, fA=self.fragmentA, fB=self.fragmentB, substructsA=ref_subnode_mols, substructsB=ref_synthon_mols,
                target=self.target, holo_file=holo_file, interactions_file=ints_file, fragment_ints_data=self.fragment_ints_data,
                mol_file=mol_file, minimised_json_file=json_file
            ) for mol, ref_subnode_mols, ref_synthon_mols, holo_file, ints_file, mol_file, json_file in
            zip(self.filter_mols, self.ref_subnode_mols, self.ref_synthon_mols, self.holo_files, self.interaction_files, self.mol_files, self.json_files)
        )
        print('Finished scoring')
        self.process_scoring_results()

    def process_scoring_results(self):
        """

        :return:
        """
        print('Processing scoring results')
        substructure_overlap_checks, energy_ratio_checks, energy_ratios, sucoses, sucos_checks = [], [], [], [], []

        passing_mols, passing_dicts, passing_sucoses = [], [], []

        for res, mol in zip(self.scoring_results, self.filter_mols):
            passing_res, data_dict = res[0], res[1]
            substructure_overlap_checks.append(data_dict['substructure_overlap_check'])
            energy_ratio_checks.append(data_dict['energy_ratio_check'])
            energy_ratios.append(data_dict['const_energy/unconst_energy'])
            sucoses.append(data_dict['sucos'])
            sucos_checks.append(data_dict['sucos_check'])

            if passing_res:
                passing_mols.append(mol)
                passing_dicts.append(data_dict)
                passing_sucoses.append(data_dict['sucos'])

        self.scoring_data.update({'substructure_overlap_check': substructure_overlap_checks,
                                  'energy_ratio_check': energy_ratio_checks,
                                  'const_energy/unconst_energy': energy_ratios,
                                  'sucos': sucoses,
                                  'sucos_check': sucos_checks})

        dump_json(self.scoring_data, os.path.join(self.scoring_dir, 'scored_data.json'))

        passing_mols = order_by_lst(passing_mols, passing_sucoses)
        passing_dicts = order_by_lst(passing_dicts, passing_sucoses)

        w = Chem.SDWriter(os.path.join(self.scoring_dir, 'ordered_mols.sdf'))
        for mol, passing_dict in zip(passing_mols, passing_dicts):
            new_mol = add_props_from_dict(mol, passing_dict)
            w.write(new_mol)
        w.close()

        print('Scoring results processed')
