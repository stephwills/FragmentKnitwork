import os
import shutil
from argparse import ArgumentParser
import time

from FragmentKnitwork.Quilter.align import alignment, placement, reverse_alignment
from FragmentKnitwork.utils import quilterConfig as config
from FragmentKnitwork.utils.quilterUtils import (get_ref_substructures,
                                                 move_files, create_alignment_dirs)
from FragmentKnitwork.utils.utils import dump_json, load_json, order_by_lst, get_mol
from joblib import Parallel, delayed
from rdkit import Chem


def main():
    parser = ArgumentParser()
    parser.add_argument('--json_file', required=False)
    parser.add_argument('--pdb_file')
    parser.add_argument('--fragmentA')
    parser.add_argument('--fragmentB')
    parser.add_argument('--substructure_dir', required=False, default=config.SUBSTRUCTURE_DIR)
    parser.add_argument('--output_dir', required=False, default=config.OUTPUT_DIR)
    parser.add_argument('--working_dir', required=False, default=config.WORKING_DIR)
    parser.add_argument('--n_cpus', type=int, required=False, default=config.N_CPUS)
    parser.add_argument('--target', required=False)
    parser.add_argument('--type_merge', choices=['pure_merge', 'impure_merge', 'reverse_merge'])
    parser.add_argument('--parallel', action='store_true', default=config.PARALLEL)
    parser.add_argument('--min_files', action='store_true', default=config.MIN_FILES)
    parser.add_argument('--move_files', action='store_true', default=config.MOVE_FILES)
    parser.add_argument('--limit_num_run', action='store_true', default=config.LIMIT_NUM_RUN)
    parser.add_argument('--max_num_filter', type=int, required=False, default=config.MAX_NUM_FILTER)
    args = parser.parse_args()

    # create pair dirs
    pair, pair_working_dir, pair_output_dir = create_alignment_dirs(args.fragmentA, args.fragmentB, args.working_dir, args.output_dir)
    if os.path.exists(os.path.join(pair_output_dir, 'passing_data.json')):
        print(f'Pair {pair} already run')
        return None

    # read in input data
    file = args.json_file
    data = load_json(file)
    smiles, cmpd_ids = [], []
    ref_subnodes, ref_synthons = [], []
    ref_subnode_mols, ref_synthon_mols = [], []

    # only used for impure merges & reverse merges
    similarities, used_subnodes, used_synthons = [], [], []

    # only used for reverse merges
    names, original_smiles, original_names, original_similarities = [], [], [], []
    rev_synthons, intermed_nodes1, intermed_nodes2 = [], [], []

    # read data from input file
    for sub_pair in data:
        print('Substructure pair', sub_pair)
        n_expans = len(data[sub_pair]['expansions'])
        print('N expansions', n_expans)

        # shared between all merge types
        smiles.extend(data[sub_pair]['expansions'])
        cmpd_ids.extend(data[sub_pair]['cmpd_ids'])

        # get ref substructures (not for reverse merge)
        subnode, synthon = sub_pair.split('*')[0], sub_pair.split('*')[1]
        ref_subnodes.extend([subnode] * n_expans)
        ref_synthons.extend([synthon] * n_expans)

        # get ref subnodes (different for reverse merge)
        if args.type_merge != 'reverse_merge':
            ref_subnode = get_ref_substructures(subnode, args.fragmentA, substructure_dir=args.substructure_dir, isSynthon=False)
            ref_subnode_mols.extend([ref_subnode] * n_expans)
        if args.type_merge == 'reverse_merge':
            ref_subnode = get_ref_substructures(subnode, args.fragmentA, substructure_dir=args.substructure_dir, isSynthon=True)
            ref_subnode_mols.extend([ref_subnode] * n_expans)

        # get ref synthons
        ref_synthon = get_ref_substructures(synthon, args.fragmentB, substructure_dir=args.substructure_dir, isSynthon=True)
        ref_synthon_mols.extend([ref_synthon] * n_expans)

        # if impure merges
        if args.type_merge == 'impure_merge' or args.type_merge == 'reverse_merge':
            similarities.extend(data[sub_pair]['similarities'])
            used_subnodes.extend([subnode] * n_expans)
            used_synthons.extend(data[sub_pair]['used_synthons'])

        if args.type_merge == 'reverse_merge':
            original_smi = sub_pair.split('*')[2]
            original_smiles.extend([original_smi] * n_expans)
            original_similarities.extend(data[sub_pair]['original_similarities'])
            original_names.extend(data[sub_pair]['original_names'])
            rev_synthons.extend(data[sub_pair]['rev_synthons'])
            intermed_nodes1.extend(data[sub_pair]['intermed_nodes1'])
            if 'intermed_nodes2' in list(data[sub_pair].keys()):
                intermed_nodes2.extend(data[sub_pair]['intermed_nodes2'])
            names.extend(data[sub_pair]['names'])

    # load mols from SMILES
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]

    # create names
    N = len(smiles)

    if args.type_merge != 'reverse_merge':
        names = [f"{args.fragmentA}-{args.fragmentB}_{idx}" for idx in range(N)]

    # limit the number of molecules to run
    if args.limit_num_run:
        N = args.max_num_filter
        if N > len(smiles):
            N = len(smiles)

    if args.type_merge == 'pure_merge':
        mols, smiles, cmpd_ids, names = mols[:N], smiles[:N], cmpd_ids[:N], names[:N]
        ref_subnodes, ref_synthons = ref_subnodes[:N], ref_synthons[:N]
        ref_subnode_mols, ref_synthon_mols = ref_subnode_mols[:N], ref_synthon_mols[:N]

    # ordered according to the similarity value of the substructue
    if args.type_merge == 'impure_merge' or args.type_merge == 'reverse_merge':
        mols = order_by_lst(mols, similarities)[:N]
        smiles = order_by_lst(smiles, similarities)[:N]
        cmpd_ids = order_by_lst(cmpd_ids, similarities)[:N]
        ref_subnodes = order_by_lst(ref_subnodes, similarities)[:N]
        ref_synthons = order_by_lst(ref_synthons, similarities)[:N]
        ref_subnode_mols = order_by_lst(ref_subnode_mols, similarities)[:N]
        ref_synthon_mols = order_by_lst(ref_synthon_mols, similarities)[:N]
        used_subnodes = order_by_lst(used_subnodes, similarities)[:N]
        used_synthons = order_by_lst(used_synthons, similarities)[:N]

        if args.type_merge == 'reverse_merge':
            original_names = order_by_lst(original_names, similarities)[:N]
            original_smiles = order_by_lst(original_smiles, similarities)[:N]
            original_similarities = order_by_lst(original_similarities, similarities)[:N]
            rev_synthons = order_by_lst(rev_synthons, similarities)[:N]
            intermed_nodes1 = order_by_lst(intermed_nodes1, similarities)[:N]
            if len(intermed_nodes2) > 0:
                intermed_nodes2 = order_by_lst(intermed_nodes2, similarities)[:N]

        similarities.sort(reverse=True)
        similarities = similarities[:N]

    data = {'names': names,
            'smiles': smiles,
            'cmpd_ids': cmpd_ids,
            'ref_subnodes': ref_subnodes,
            'ref_synthons': ref_synthons}

    if args.type_merge == 'impure_merge' or args.type_merge == 'reverse_merge':
        data.update({'similarities': similarities,
                     'used_subnodes': used_subnodes,
                     'used_synthons': used_synthons})

    if args.type_merge == 'reverse_merge':
        data.update({'original_smiles': original_smiles,
                     'original_names': original_names,
                     'original_similarities': original_similarities,
                     'rev_synthons': rev_synthons,
                     'intermed_nodes1': intermed_nodes1})
        if len(intermed_nodes2) > 0:
            data['intermed_nodes2'] = intermed_nodes2

    # write the input data
    dump_json(data, os.path.join(pair_output_dir, 'input_data.json'))
    print('Input data dumped')

    # run the filtering
    start = time.time()
    if args.type_merge == 'impure_merge':
        results = Parallel(n_jobs=args.n_cpus, backend='multiprocessing')(
                delayed(alignment)(
                    smi, name, merge, ref_subnode, ref_synthon, used_subnode, used_synthon,
                    args.pdb_file, pair_working_dir
                ) for smi, name, merge, ref_subnode, ref_synthon, used_subnode, used_synthon in
                zip(smiles, names, mols, ref_subnode_mols, ref_synthon_mols, used_subnodes, used_synthons)
            )
    if args.type_merge == 'pure_merge':
        results = Parallel(n_jobs=args.n_cpus, backend='multiprocessing')(
                delayed(placement)(
                    smi, name, ref_subnode, ref_synthon, args.pdb_file, pair_working_dir
                ) for smi, name, ref_subnode, ref_synthon, in
                zip(smiles, names, ref_subnode_mols, ref_synthon_mols)
            )
    if args.type_merge == 'reverse_merge':
        if len(intermed_nodes2) == 0:
            results = Parallel(n_jobs=args.n_cpus, backend='multiprocessing')(
                    delayed(reverse_alignment)(
                        smi, original_smi, intermed_smi1, name, merge, ref_subnode, ref_synthon, rev_synthon, used_synthon,
                        args.pdb_file, pair_working_dir
                    ) for smi, original_smi, intermed_smi1, name, merge, ref_subnode, ref_synthon, rev_synthon, used_synthon in
                    zip(smiles, original_smiles, intermed_nodes1, names, mols, ref_subnode_mols, ref_synthon_mols, rev_synthons, used_synthons)
                )
        else:
            results = Parallel(n_jobs=args.n_cpus, backend='multiprocessing')(
                    delayed(reverse_alignment)(
                        smi, original_smi, intermed_smi1, name, merge, ref_subnode, ref_synthon, rev_synthon, used_synthon,
                        args.pdb_file, pair_working_dir, intermed_smi2
                    ) for smi, original_smi, intermed_smi1, name, merge, ref_subnode, ref_synthon, rev_synthon, used_synthon, intermed_smi2 in
                    zip(smiles, original_smiles, intermed_nodes1, names, mols, ref_subnode_mols, ref_synthon_mols, rev_synthons, used_synthons, intermed_nodes2)
                )
    end = time.time()
    total_time = round(end-start, 2)

    # record the output data
    passes = [res[0] for res in results]
    aligned_mols = [res[1] for res in results]
    results_dirs = [res[2] for res in results]
    errors = [res[3] for res in results]
    timings = [res[4] for res in results]
    data['passes'] = passes
    data['errors'] = errors
    data['timings'] = timings

    # record overall timing data
    data['total_time'] = total_time
    data['n_cpus'] = args.n_cpus

    # write the output data
    dump_json(data, os.path.join(pair_output_dir, 'passing_data.json'))

    # write an SDF with the aligned molecules
    w = Chem.SDWriter(os.path.join(pair_output_dir, f"aligned_mols.sdf"))
    n_aligned = 0
    if args.type_merge == 'pure_merge':
        similarities = [None] * len(names)
    for name, smi, cmpd_id, similarity, mol in zip(names, smiles, cmpd_ids, similarities, aligned_mols):
        if mol:
            mol.SetProp('_Name', name)
            mol.SetProp('smiles', smi)
            mol.SetProp('cmpd_ids', ','.join(cmpd_id))
            if args.type_merge != 'pure_merge':
                mol.SetProp('similarity', str(similarity))
            w.write(mol)
            n_aligned += 1
    print(n_aligned, 'conformers generated out of', len(smiles))

    # move files out of the work directory
    if args.move_files:
        results_dirs = [i for i in results_dirs if i]
        move_files(results_dirs, pair_output_dir, args.min_files)
        shutil.rmtree(pair_working_dir)


if __name__ == "__main__":
    main()
