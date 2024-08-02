import os
import time
from argparse import ArgumentParser

from FragmentKnitwork.Quilter.align import r_group_alignment
from FragmentKnitwork.utils import quilterConfig as config
from FragmentKnitwork.utils.utils import dump_json, get_protein, get_prop_as_list
from joblib import Parallel, delayed
from rdkit import Chem


def main():
    """
    Run R-group alignment using an sdf of placed merges

    :return:
    """
    parser = ArgumentParser()
    parser.add_argument('--sdf_file', required=True)
    parser.add_argument('--original_output_dir', required=True)
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--n_cpus', type=int, required=False, default=config.N_CPUS)
    parser.add_argument('--target', required=True)
    parser.add_argument('--type_merge')
    parser.add_argument('--parallel', action='store_true', default=config.PARALLEL)
    parser.add_argument('--min_files', action='store_true', default=config.MIN_FILES)
    parser.add_argument('--move_files', action='store_true', default=config.MOVE_FILES)
    parser.add_argument('--limit_num_run', action='store_true', default=config.LIMIT_NUM_RUN)
    parser.add_argument('--max_num_filter', type=int, required=False, default=config.MAX_NUM_FILTER)
    args = parser.parse_args()

    mols = list(Chem.SDMolSupplier(args.sdf_file))

    # read in data using the mol properties
    names = get_prop_as_list(mols, '_Name')
    smiles = get_prop_as_list(mols, 'smiles')
    fAs, fBs = get_prop_as_list(mols, 'fragmentA'), get_prop_as_list(mols, 'fragmentB')
    pdb_files = [get_protein(args.target, fA, protonated=False, desolv=True) for fA in fAs]
    similarity = get_prop_as_list(mols, 'similarity')
    og_cmpd_id = get_prop_as_list(mols, 'cmpd_id')
    used_subnodes, used_synthons = get_prop_as_list(mols, 'used_subnode'), get_prop_as_list(mols, 'used_synthon')
    ref_subnodes, ref_synthons = get_prop_as_list(mols, 'ref_subnode'), get_prop_as_list(mols, 'ref_synthon')
    r_groups = get_prop_as_list(mols, 'r_group')
    smiles_before_rgrp_expansions = get_prop_as_list(mols, 'smiles_before_rgrp_expansions')

    data = {
        'fragmentA': fAs,
        'fragmentB': fBs,
        'names': names,
        'smiles': smiles,
        'similarity': similarity,
        'og_cmpd_id': og_cmpd_id,
        'used_subnode': used_subnodes,
        'used_synthon': used_synthons,
        'ref_subnode': ref_subnodes,
        'ref_synthon': ref_synthons,
        'r_group': r_groups,
        'smiles_before_rgrp_expansions': smiles_before_rgrp_expansions
    }

    if args.type_merge == 'reverse_merge':
        data['rev_synthon'] = get_prop_as_list(mols, 'rev_synthon')
        data['intermed_node1'] = get_prop_as_list(mols, 'intermed_node1')
    dump_json(data, os.path.join(args.output_dir, 'input_data.json'))

    # get original molecules for performing alignment with
    og_mols = []
    for fA, fB, name in zip(fAs, fBs, names):
        og_name = name.split('_rgrp')[0]
        og_name_hyph = og_name.replace('_', '-')
        mol_file = os.path.join(args.original_output_dir, f"{fA}-{fB}", f"{og_name_hyph}", f"{og_name_hyph}.minimised.mol")
        og_mol = Chem.MolFromMolFile(mol_file)
        og_mols.append(og_mol)

    start = time.time()
    results = Parallel(n_jobs=args.n_cpus, backend='multiprocessing')(
        delayed(r_group_alignment)(
            smi, og_mol, args.output_dir, pdb_file, name
        ) for smi, og_mol, pdb_file, name in zip(smiles, og_mols, pdb_files, names)
    )
    end = time.time()
    time_taken = round(end-start, 2)
    print('Time taken:', time_taken)

    # write the output data
    passes = [res[0] for res in results]
    timings = [res[2] for res in results]
    data['passes'] = passes
    data['timings'] = timings
    data['time'] = time_taken
    data['n_cpus'] = args.n_cpus
    dump_json(data, os.path.join(args.output_dir, 'passing_data.json'))


if __name__ == "__main__":
    main()
