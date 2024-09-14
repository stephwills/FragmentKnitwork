import os
import shutil
import time
from FragmentKnitwork.Quilter.alignWictor import (victor_alignment,
                                                  wictor_alignment,
                                                  wictor_placement)
from FragmentKnitwork.utils import quilterConfig as config
from FragmentKnitwork.utils.quilterUtils import (create_alignment_dirs,
                                                 get_ref_substructures,
                                                 move_files, add_props_from_dict)
from FragmentKnitwork.utils.fragmentUtils import filt_by_bool
from FragmentKnitwork.utils.utils import dump_json, load_json, order_by_lst
from FragmentKnitwork.Quilter.checkPlacedMol import process_molecule
from joblib import Parallel, delayed
from rdkit import Chem


def write_ref_substructures(pair_output_dir, best_ref_subnodes, best_ref_synthons):
    """

    :param pair_output_dir:
    :param best_ref_subnodes:
    :param best_ref_synthons:
    :return:
    """
    w1 = Chem.SDWriter(os.path.join(pair_output_dir, f"wictor_ref_subnodes.sdf"))
    w2 = Chem.SDWriter(os.path.join(pair_output_dir, f"wictor_ref_synthons.sdf"))

    n_aligned = 0
    dummy_mol = Chem.MolFromSmiles('')

    for best_ref_subnode, best_ref_synthon in zip(best_ref_subnodes, best_ref_synthons):
        if best_ref_subnode and best_ref_synthon:
            w1.write(best_ref_subnode)
            w2.write(best_ref_synthon)
            n_aligned += 1
        else:
            w1.write(dummy_mol)
            w2.write(dummy_mol)
    print(n_aligned, 'alignments generated using wictor out of', len(best_ref_subnodes))


class WictorPipeline:

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
                 limit_num_minimize=None,
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
        self.limit_num_minimize = limit_num_minimize
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

        if self.merge_type == 'impure_merge':
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

        if self.merge_type == 'impure_merge':
            data.update({'similarities': self.similarities,
                        'used_subnodes': self.used_subnodes,
                        'used_synthons': self.used_synthons})

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
            ref_subnode = get_ref_substructures(subnode, self.fragmentA, substructure_dir=self.substructure_dir, isSynthon=False)
            self.ref_subnode_mols.append(ref_subnode)
            ref_synthon = get_ref_substructures(synthon, self.fragmentB, substructure_dir=self.substructure_dir, isSynthon=True)
            self.ref_synthon_mols.append(ref_synthon)

    def run_wictor_filtering(self):
        """

        :return:
        """
        start = time.time()
        if self.merge_type == 'pure_merge':
            self.wictor_results = Parallel(n_jobs=self.n_cpus, backend='multiprocessing')(
                delayed(wictor_placement)(
                    smi, name, ref_subnode, ref_synthon, self.pdb_file, self.pair_working_dir
                ) for smi, name, ref_subnode, ref_synthon, in
                zip(self.smiles, self.names, self.ref_subnode_mols, self.ref_synthon_mols)
            )
        if self.merge_type == 'impure_merge':
            self.wictor_results = Parallel(n_jobs=self.n_cpus, backend='multiprocessing')(
                delayed(wictor_alignment)(
                    smi, name, merge, ref_subnode, ref_synthon, used_subnode, used_synthon,
                    self.pdb_file, self.pair_working_dir
                ) for smi, name, merge, ref_subnode, ref_synthon, used_subnode, used_synthon in
                zip(self.smiles, self.names, self.mols, self.ref_subnode_mols, self.ref_synthon_mols, self.used_subnodes, self.used_synthons)
            )
        end = time.time()
        self.total_wictor_time = round(end - start, 2)

    def save_wictor_results(self):
        """

        :return:
        """
        # results contain: passes, sucos_scores, best_ref_subnodes, best_ref_synthons, mappings, wictor_timings
        print('Processing filter results')
        self.wictor_passes, self.wictor_sucoses, self.best_ref_subnodes, self.best_ref_synthons, self.mappings, self.wictor_timings = [], [], [], [], [], []
        for res in self.wictor_results:
            self.wictor_passes.append(res[0])
            self.wictor_sucoses.append(res[1])
            self.best_ref_subnodes.append(res[2])
            self.best_ref_synthons.append(res[3])
            self.mappings.append(res[4])
            self.wictor_timings.append(res[5])

        # write summary of results to file
        self.data.update({'passes': self.wictor_passes,
                          'wictor_sucoses': self.wictor_sucoses,
                          'wictor_timings': self.wictor_timings,
                          'total_wictor_time': self.total_wictor_time,
                          'n_cpus': self.n_cpus})

        dump_json(self.data, os.path.join(self.pair_output_dir, 'wictor_data.json'))

    def set_max_minimize(self):
        """

        :return:
        """
        self.names = order_by_lst(self.names, self.wictor_sucoses)[:self.limit_num_minimize]
        self.mols = order_by_lst(self.mols, self.wictor_sucoses)[:self.limit_num_minimize]
        self.smiles = order_by_lst(self.smiles, self.wictor_sucoses)[:self.limit_num_minimize]
        self.cmpd_ids = order_by_lst(self.cmpd_ids, self.wictor_sucoses)[:self.limit_num_minimize]
        self.ref_subnodes = order_by_lst(self.ref_subnodes, self.wictor_sucoses)[:self.limit_num_minimize]
        self.ref_synthons = order_by_lst(self.ref_synthons, self.wictor_sucoses)[:self.limit_num_minimize]
        self.ref_subnode_mols = order_by_lst(self.ref_subnode_mols, self.wictor_sucoses)[:self.limit_num_minimize]
        self.ref_synthon_mols = order_by_lst(self.ref_synthon_mols, self.wictor_sucoses)[:self.limit_num_minimize]
        self.best_ref_subnodes = order_by_lst(self.best_ref_subnodes, self.wictor_sucoses)[:self.limit_num_minimize]
        self.best_ref_synthons = order_by_lst(self.best_ref_synthons, self.wictor_sucoses)[:self.limit_num_minimize]
        self.mappings = order_by_lst(self.mappings, self.wictor_sucoses)[:self.limit_num_minimize]

        if self.merge_type == 'impure_merge':
            self.used_subnodes = order_by_lst(self.used_subnodes, self.wictor_sucoses)[:self.limit_num_minimize]
            self.used_synthons = order_by_lst(self.used_synthons, self.wictor_sucoses)[:self.limit_num_minimize]
            self.similarities = order_by_lst(self.similarities, self.wictor_sucoses)[:self.limit_num_minimize]

        self.wictor_sucoses.sort(reverse=True)
        self.wictor_sucoses = self.wictor_sucoses[:self.limit_num_minimize]

    def process_wictor_results(self):
        """

        :return:
        """
        self.save_wictor_results()

        # write an SDF with the ref subnodes and synthons
        write_ref_substructures(self.pair_output_dir, self.best_ref_subnodes, self.best_ref_synthons)

        # filter according to the filter results
        self.names = filt_by_bool(self.names, self.wictor_passes)
        self.mols = filt_by_bool(self.mols, self.wictor_passes)
        self.smiles = filt_by_bool(self.smiles, self.wictor_passes)
        self.cmpd_ids = filt_by_bool(self.cmpd_ids, self.wictor_passes)
        self.ref_subnodes = filt_by_bool(self.ref_subnodes, self.wictor_passes)
        self.ref_synthons = filt_by_bool(self.ref_synthons, self.wictor_passes)
        self.ref_subnode_mols = filt_by_bool(self.ref_subnode_mols, self.wictor_passes)
        self.ref_synthon_mols = filt_by_bool(self.ref_synthon_mols, self.wictor_passes)

        if self.merge_type == 'impure_merge':
            self.used_subnodes = filt_by_bool(self.used_subnodes, self.wictor_passes)
            self.used_synthons = filt_by_bool(self.used_synthons, self.wictor_passes)
            self.similarities = filt_by_bool(self.similarities, self.wictor_passes)

        self.best_ref_subnodes = filt_by_bool(self.best_ref_subnodes, self.wictor_passes)
        self.best_ref_synthons = filt_by_bool(self.best_ref_synthons, self.wictor_passes)
        self.mappings = filt_by_bool(self.mappings, self.wictor_passes)
        self.wictor_sucoses = filt_by_bool(self.wictor_sucoses, self.wictor_passes)

        if self.limit_num_minimize:
            self.set_max_minimize()

        self.wictor_data = self.get_data_dict()

    def run_filtering(self):
        """
        Run the relevant filtering/conformer generation method depending on the type of merge

        :return:
        """
        self.run_wictor_filtering()
        self.process_wictor_results()
        self.run_victor_filtering()
        self.process_victor_results()

    def run_victor_filtering(self):
        """

        :return:
        """
        start = time.time()
        self.victor_results = Parallel(n_jobs=self.n_cpus, backend='multiprocessing')(
            delayed(victor_alignment)(
                smi, name, ref_subnode, ref_synthon, mapping, self.pdb_file, self.pair_working_dir
            ) for smi, name, ref_subnode, ref_synthon, mapping in
            zip(self.smiles,self.names, self.best_ref_subnodes, self.best_ref_synthons, self.mappings)
        )
        end = time.time()
        self.total_victor_time = round(end - start, 2)

    def save_victor_results(self):
        """

        :return:
        """
        # results contain: passes, aligned_mols, results_dirs, errors, timings
        print('Processing filter results')
        self.filter_passes, self.filter_mols, self.results_dirs, self.errors, self.timings = [], [], [], [], []
        for res in self.victor_results:
            self.filter_passes.append(res[0])
            self.filter_mols.append(res[1])
            self.results_dirs.append(res[2])
            self.errors.append(res[3])
            self.timings.append(res[4])

        # write summary of results to file
        self.wictor_data.update({'passes': self.filter_passes,
                                 'errors': self.errors,
                                 'timings': self.timings,
                                 'total_filter_time': self.total_victor_time,
                                 'n_cpus': self.n_cpus})

        dump_json(self.wictor_data, os.path.join(self.pair_output_dir, 'passing_data.json'))

        print('Moving successful results files to output dir')
        self.move_results_files()

    def process_victor_results(self):
        """

        :return:
        """
        # save the results to output json file (and basic processing)
        self.save_victor_results()

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

        if self.merge_type == 'impure_merge':
            self.used_subnodes = filt_by_bool(self.used_subnodes, self.filter_passes)
            self.used_synthons = filt_by_bool(self.used_synthons, self.filter_passes)
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

        dump_json(self.scoring_data, os.path.join(self.scoring_dir, f'{self.fragmentA}-{self.fragmentB}_scored_data.json'))

        passing_mols = order_by_lst(passing_mols, passing_sucoses)
        passing_dicts = order_by_lst(passing_dicts, passing_sucoses)

        w = Chem.SDWriter(os.path.join(self.scoring_dir, f'{self.fragmentA}-{self.fragmentB}_ordered_mols.sdf'))
        for mol, passing_dict in zip(passing_mols, passing_dicts):
            new_mol = add_props_from_dict(mol, passing_dict)
            w.write(new_mol)
        w.close()

        print('Scoring results processed')
