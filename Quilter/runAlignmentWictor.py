import os
import shutil
import time
from argparse import ArgumentParser

from FragmentKnitwork.Quilter.alignWictor import (victor_alignment,
                                                  wictor_alignment,
                                                  wictor_placement)
from FragmentKnitwork.Quilter.scoreAlignedMols import process_molecule
from FragmentKnitwork.utils import quilterConfig as config
from FragmentKnitwork.utils.quilterUtils import (add_props_from_dict,
                                                 create_alignment_dirs,
                                                 get_ref_substructures,
                                                 move_files)
from FragmentKnitwork.utils.utils import dump_json, load_json, order_by_lst
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


def get_data_for_minimization(passes, smiles, names, cmpd_ids, sucos_scores, best_ref_subnodes, best_ref_synthons,
                              mappings, max_num_minimize=None):
    """

    :param passes:
    :param smiles:
    :param names:
    :param cmpd_ids:
    :param sucos_scores:
    :param best_ref_subnodes:
    :param best_ref_synthons:
    :param mappings:
    :param max_num_minimize:
    :return:
    """
    minimize_smiles, minimize_names, minimize_ref_subnodes, minimize_ref_synthons = [], [], [], []
    minimize_mappings, minimize_scores, minimize_cmpd_ids = [], [], []

    for res, smi, name, cmpd_id, score, ref_subnode, ref_synthon, mapping in zip(passes, smiles, names, cmpd_ids,
                                                                                 sucos_scores, best_ref_subnodes,
                                                                                 best_ref_synthons, mappings):
        if res:
            minimize_smiles.append(smi)
            minimize_names.append(name)
            minimize_ref_subnodes.append(ref_subnode)
            minimize_ref_synthons.append(ref_synthon)
            minimize_mappings.append(mapping)
            minimize_scores.append(score)
            minimize_cmpd_ids.append(cmpd_id)

    if max_num_minimize:
        Nmin = max_num_minimize
        if Nmin > len(minimize_smiles):
            Nmin = len(minimize_smiles)

    else:
        Nmin = len(minimize_smiles)

    # order by sucos, select top X
    minimize_smiles = order_by_lst(minimize_smiles, minimize_scores)[:Nmin]
    minimize_names = order_by_lst(minimize_names, minimize_scores)[:Nmin]
    minimize_cmpd_ids = order_by_lst(minimize_cmpd_ids, minimize_scores)[:Nmin]
    minimize_ref_subnodes = order_by_lst(minimize_ref_subnodes, minimize_scores)[:Nmin]
    minimize_ref_synthons = order_by_lst(minimize_ref_synthons, minimize_scores)[:Nmin]
    minimize_mappings = order_by_lst(minimize_mappings, minimize_scores)[:Nmin]
    minimize_scores.sort(reverse=True)
    minimize_scores = minimize_scores[:Nmin]

    return minimize_smiles, minimize_names, minimize_cmpd_ids, minimize_ref_subnodes, minimize_ref_synthons, minimize_mappings, minimize_scores


def add_wictor_dict_info_to_victor_dict(data, minimize_data):
    """

    :param data:
    :param minimize_data:
    :return:
    """
    # add other data from wictor_data
    wictor_names = data['names']
    minimize_names = minimize_data['names']

    # wictor lists not in victor
    wictor_keys = list(data.keys())
    minimize_keys = list(minimize_data.keys())
    keys_get = [i for i in wictor_keys if i not in minimize_keys]

    # get idxs of mols that have been maintained
    idxs = [wictor_names.index(name) for name in minimize_names]

    for k in keys_get:
        if 'time' not in k:
            new_lst = [data[k][idx] for idx in idxs]
            minimize_data[k] = new_lst

    return minimize_data


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
    parser.add_argument('--type_merge', choices=['pure_merge', 'impure_merge'])
    parser.add_argument('--parallel', action='store_true', default=config.PARALLEL)
    parser.add_argument('--min_files', action='store_true', default=config.MIN_FILES)
    parser.add_argument('--move_files', action='store_true', default=config.MOVE_FILES)
    parser.add_argument('--limit_num_run', action='store_true')
    parser.add_argument('--max_num_filter', type=int, required=False, default=None)
    parser.add_argument('--max_num_minimize', type=int, required=False, default=None)
    parser.add_argument('--run_scoring', action='store_true')
    args = parser.parse_args()
    start = time.time()

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

    # read data from file
    for sub_pair in data:
        print('Substructure pair', sub_pair)
        n_expans = len(data[sub_pair]['expansions'])
        print('N expansions', n_expans)

        # shared between all merge types
        smiles.extend(data[sub_pair]['expansions'])
        cmpd_ids.extend(data[sub_pair]['cmpd_ids'])

        # get ref substructures
        subnode, synthon = sub_pair.split('*')[0], sub_pair.split('*')[1]
        ref_subnodes.extend([subnode] * n_expans)
        ref_synthons.extend([synthon] * n_expans)

        # get the ref subnode and synthon molecules
        ref_subnode = get_ref_substructures(subnode, args.fragmentA, substructure_dir=args.substructure_dir, isSynthon=False)
        ref_subnode_mols.extend([ref_subnode] * n_expans)
        ref_synthon = get_ref_substructures(synthon, args.fragmentB, substructure_dir=args.substructure_dir, isSynthon=True)
        ref_synthon_mols.extend([ref_synthon] * n_expans)

        # if impure merges
        if args.type_merge == 'impure_merge':
            similarities.extend(data[sub_pair]['similarities'])
            used_subnodes.extend([subnode] * n_expans)
            used_synthons.extend(data[sub_pair]['used_synthons'])

    # load mols from SMILES
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]

    # create names
    N = len(smiles)
    names = [f"{args.fragmentA}-{args.fragmentB}_{idx}" for idx in range(N)]

    if args.limit_num_run:
        N = args.max_num_filter
        if N > len(smiles):
            N = len(smiles)

    if args.type_merge == 'pure_merge':
        mols, smiles, cmpd_ids, names = mols[:N], smiles[:N], cmpd_ids[:N], names[:N]
        ref_subnodes, ref_synthons = ref_subnodes[:N], ref_synthons[:N]
        ref_subnode_mols, ref_synthon_mols = ref_subnode_mols[:N], ref_synthon_mols[:N]

    if args.type_merge == 'impure_merge':
        mols = order_by_lst(mols, similarities)[:N]
        smiles = order_by_lst(smiles, similarities)[:N]
        cmpd_ids = order_by_lst(cmpd_ids, similarities)[:N]
        ref_subnodes = order_by_lst(ref_subnodes, similarities)[:N]
        ref_synthons = order_by_lst(ref_synthons, similarities)[:N]
        ref_subnode_mols = order_by_lst(ref_subnode_mols, similarities)[:N]
        ref_synthon_mols = order_by_lst(ref_synthon_mols, similarities)[:N]
        used_subnodes = order_by_lst(used_subnodes, similarities)[:N]
        used_synthons = order_by_lst(used_synthons, similarities)[:N]

        similarities.sort(reverse=True)
        similarities = similarities[:N]

    data = {'names': names,
            'smiles': smiles,
            'cmpd_ids': cmpd_ids,
            'ref_subnodes': ref_subnodes,
            'ref_synthons': ref_synthons}

    if args.type_merge == 'impure_merge':
        data.update({'similarities': similarities,
                     'used_subnodes': used_subnodes,
                     'used_synthons': used_synthons})

    # write the input data
    dump_json(data, os.path.join(pair_output_dir, 'input_data.json'))
    print('Input data dumped')

    ############################## run alignment with wictor ##############################
    # run the filtering
    start = time.time()
    if args.type_merge == 'impure_merge':
        results = Parallel(n_jobs=args.n_cpus, backend='multiprocessing')(
                delayed(wictor_alignment)(
                    smi, name, merge, ref_subnode, ref_synthon, used_subnode, used_synthon,
                    args.pdb_file, pair_working_dir
                ) for smi, name, merge, ref_subnode, ref_synthon, used_subnode, used_synthon in
                zip(smiles, names, mols, ref_subnode_mols, ref_synthon_mols, used_subnodes, used_synthons)
            )
    if args.type_merge == 'pure_merge':
        results = Parallel(n_jobs=args.n_cpus, backend='multiprocessing')(
                delayed(wictor_placement)(
                    smi, name, ref_subnode, ref_synthon, args.pdb_file, pair_working_dir
                ) for smi, name, ref_subnode, ref_synthon, in
                zip(smiles, names, ref_subnode_mols, ref_synthon_mols)
            )
    end = time.time()
    total_time = round(end-start, 2)

    # record the output data
    passes = [res[0] for res in results]
    sucos_scores = [res[1] for res in results]
    best_ref_subnodes = [res[2] for res in results]
    best_ref_synthons = [res[3] for res in results]
    mappings = [res[4] for res in results]
    wictor_timings = [res[5] for res in results]
    data['passes'] = passes
    data['sucoses'] = sucos_scores
    data['mappings'] = mappings
    data['total_wictor_time'] = total_time
    data['wictor_timings'] = wictor_timings
    data['n_cpus'] = args.n_cpus

    # dumpthe wictor output data
    wictor_data_fname = os.path.join(pair_output_dir, 'wictor_data.json')
    dump_json(data, wictor_data_fname)

    # write an SDF with the ref subnodes and synthons
    write_ref_substructures(pair_output_dir, best_ref_subnodes, best_ref_synthons)

    #######################################################################################

    minimize_smiles, minimize_names, minimize_cmpd_ids, minimize_ref_subnodes, minimize_ref_synthons, minimize_mappings, minimize_scores = get_data_for_minimization(
        passes, smiles, names, cmpd_ids, sucos_scores, best_ref_subnodes, best_ref_synthons, mappings, max_num_minimize=args.max_num_minimize)

    minimize_data = {'smiles': minimize_smiles,
                     'names': minimize_names,
                     'cmpd_ids': minimize_cmpd_ids,
                     'mappings': minimize_mappings,
                     'wictor_sucos': minimize_scores}

    # record minimize time
    minimize_start = time.time()

    ############################## run alignment with victor ##############################

    minimize_results = Parallel(n_jobs=args.n_cpus, backend='multiprocessing')(
        delayed(victor_alignment)(
            smi, name, ref_subnode, ref_synthon, mapping, args.pdb_file, pair_working_dir
        ) for smi, name, ref_subnode, ref_synthon, mapping in
        zip(minimize_smiles, minimize_names, minimize_ref_subnodes, minimize_ref_synthons, minimize_mappings)
    )

    # record minimize time
    minimize_end = time.time()
    total_minimize = round(minimize_end - minimize_start, 2)

    # write output data for the minimization
    minimize_passes = [res[0] for res in minimize_results]
    aligned_mols = [res[1] for res in minimize_results]
    results_dirs = [res[2] for res in minimize_results]
    errors = [res[3] for res in minimize_results]
    victor_timings = [res[4] for res in minimize_results]

    minimize_data['passes'] = minimize_passes
    minimize_data['errors'] = errors
    minimize_data['total_minimize_time'] = total_minimize
    minimize_data['victor_timings'] = victor_timings
    minimize_data['n_cpus'] = args.n_cpus

    # get the other relevant data from wictor data
    minimize_data = add_wictor_dict_info_to_victor_dict(data, minimize_data)

    # write new output file for victor placements
    dump_json(minimize_data, os.path.join(pair_output_dir, 'passing_data.json'))

    # move files out of the work directory
    if args.move_files:
        results_dirs = [i for i in results_dirs if i]
        move_files(results_dirs, pair_output_dir, args.min_files)
        shutil.rmtree(pair_working_dir)

    # write an SDF with the aligned molecules
    w = Chem.SDWriter(os.path.join(pair_output_dir, f"aligned_mols.sdf"))
    n_aligned = 0
    for name, smi, cmpd_id, mol in zip(minimize_names, minimize_smiles, minimize_cmpd_ids, aligned_mols):
        if mol:
            mol.SetProp('_Name', name)
            mol.SetProp('smiles', smi)
            mol.SetProp('cmpd_ids', ','.join(cmpd_id))
            w.write(mol)
            n_aligned += 1
    print(n_aligned, 'conformers generated out of total', len(smiles), 'and out of', len(minimize_smiles), 'total minimized')

    end = time.time()
    print('Full time taken:', round(end-start, 2), 'seconds')

    #######################################################################################

    if args.run_scoring:
        print('Running scoring')
        passing_names_nohyph = [name for name, res in zip(minimize_names, passes) if res]
        passing_names = [name.replace('_', '-') for name in passing_names_nohyph]
        passing_smiles = [smi for smi, res in zip(minimize_smiles, passes) if res]
        passing_cmpd_ids = [cmpd_id for cmpd_id, res in zip(minimize_cmpd_ids, passes) if res]
        passing_mols = [mol for mol in aligned_mols if mol]

        holo_files = [os.path.join(pair_output_dir, name, f"{name}.holo_minimised.pdb") for name in passing_names]
        mol_files = [os.path.join(pair_output_dir, name, f"{name}.minimised.mol") for name in passing_names]
        json_files = [os.path.join(pair_output_dir, name, f"{name}.minimised.json") for name in passing_names]
        interactions_files = [os.path.join(pair_output_dir, name, f"{name}_interactions.json") for name in passing_names]

        fragment_ints_file = os.path.join(args.substructure_dir, 'fragment_interactions.json')
        fragment_ints_data = load_json(fragment_ints_file)

        scoring_checks = Parallel(n_jobs=args.n_cpus, backend="multiprocessing")(
            delayed(process_molecule)(
                mol, args.fragmentA, args.fragmentB, [ref_subnode], [ref_synthon], args.target, holo_file, int_file, fragment_ints_data, mol_file, json_file
            ) for mol, ref_subnode, ref_synthon, holo_file, int_file, mol_file, json_file in zip(passing_mols, minimize_ref_subnodes, minimize_ref_synthons,
                                                                                                 holo_files, interactions_files, mol_files, json_files)
        )

        wM = Chem.SDWriter(os.path.join(pair_output_dir, 'filtered_mols.sdf'))
        wA = Chem.SDWriter(os.path.join(pair_output_dir, 'filtered_fA_mols.sdf'))
        wB = Chem.SDWriter(os.path.join(pair_output_dir, 'filtered_fB_mols.sdf'))

        for name, smi, cmpd_id, mol, check in zip(passing_names_nohyph, passing_smiles, passing_cmpd_ids, scoring_checks):
            filter, mol_info, fA_mol, fB_mol = check[0], check[1], check[2], check[3]
            mol = add_props_from_dict(mol, mol_info)
            mol.SetProp('_Name', name)
            mol.SetProp('smiles', smi)
            mol.SetProp('cmpd_ids', ','.join(cmpd_id))
            wM.write(mol)
            wA.write(fA_mol)
            wB.write(fB_mol)

        wM.close()
        wA.close()
        wB.close()


if __name__ == "__main__":
    main()
