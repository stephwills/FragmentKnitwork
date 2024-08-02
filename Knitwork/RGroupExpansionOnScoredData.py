"""Perform R group expansion of a given merge"""

import os
from argparse import ArgumentParser

from FragmentKnitwork.Knitwork.queries import single_expansion
from FragmentKnitwork.utils.utils import get_intersect, load_json, dump_json
from FragmentKnitwork.utils.quilterUtils import add_props_to_new_mol
from joblib import Parallel, delayed
from rdkit import Chem
from tqdm import tqdm
import time


def get_r_group_expansions(merge, _ref_subnode, ref_synthon, fragmentA, fragmentB, r_group_data, equiv_synthon_data):
    if '[Xe]' in _ref_subnode:
        ref_subnodes = [_ref_subnode]
    else:
        if _ref_subnode in equiv_synthon_data[fragmentA]:
            ref_subnodes = equiv_synthon_data[fragmentA][_ref_subnode].split(',')
        else:
            print(_ref_subnode, 'not in equiv synthon dict')
            return [], [], []
    fA_expansions = []
    fB_expansions = []
    for ref_subnode in ref_subnodes:
        if ref_subnode in r_group_data[fragmentA]:
            _fA_expansions = r_group_data[fragmentA][ref_subnode]
            fA_expansions.extend(_fA_expansions)
    if ref_synthon in r_group_data[fragmentB]:
        fB_expansions = r_group_data[fragmentB][ref_synthon]

    r_groups = fA_expansions + fB_expansions
    r_groups = list(set(r_groups))

    all_expansions = []
    all_cmpd_ids = []
    all_r_groups = []

    if len(r_groups) > 0:
        for r_group in r_groups:
            expansions, cmpd_ids = single_expansion(merge, r_group)
            if len(expansions) > 0:
                all_expansions.extend(expansions)
                all_cmpd_ids.extend(cmpd_ids)
                all_r_groups.extend([r_group] * len(expansions))

    if len(all_expansions) > 0:
        print(f"{merge}: {len(all_expansions)} R group expansions found")
    return all_expansions, all_cmpd_ids, all_r_groups


def main():
    parser = ArgumentParser()
    parser.add_argument('--sdf_file')
    parser.add_argument('--equivalent_subnode_file')
    parser.add_argument('--r_group_expansion_file')
    parser.add_argument('--merge_type')
    parser.add_argument('--name_tag')
    parser.add_argument('--target')
    parser.add_argument('--n_cpus', type=int, required=False)
    parser.add_argument('--parallel', action='store_true')
    parser.add_argument('--output_dir')
    args = parser.parse_args()

    mols = list(Chem.SDMolSupplier(args.sdf_file))
    smiles = [mol.GetProp('smiles') for mol in mols]
    fAs = [mol.GetProp('fragmentA') for mol in mols]
    fBs = [mol.GetProp('fragmentB') for mol in mols]
    ref_subnodes = [mol.GetProp('ref_subnode') for mol in mols]
    ref_synthons = [mol.GetProp('ref_synthon') for mol in mols]

    r_group_data = load_json(args.r_group_expansion_file)
    equiv_synthon_data = load_json(args.equivalent_subnode_file)

    timings_file = os.path.join(args.output_dir, 'timings.json')
    start = time.time()
    if not args.parallel:
        for smi, fA, fB, ref_subnode, ref_synthon in tqdm(zip(smiles, fAs, fBs, ref_subnodes, ref_synthons)):
            results = get_r_group_expansions(smi, ref_subnode, ref_synthon, fA, fB, r_group_data, equiv_synthon_data)

    else:
        results = Parallel(n_jobs=args.n_cpus, backend='multiprocessing')(
            delayed(get_r_group_expansions)(
                smi, ref_subnode, ref_synthon, fA, fB,r_group_data, equiv_synthon_data
            ) for smi, fA, fB, ref_subnode, ref_synthon in
            tqdm(zip(smiles, fAs, fBs, ref_subnodes, ref_synthons))
        )

    print('Querying done, starting processing')
    new_mols = []
    for mol, res in zip(mols, results):
        all_expansions, all_cmpd_ids, all_r_groups = res[0], res[1], res[2]
        if len(all_expansions) > 0:
            for i, (expansion, _cmpd_id, r_group) in enumerate(zip(all_expansions, all_cmpd_ids, all_r_groups)):
                cmpd_id = ','.join(_cmpd_id)
                old_name = mol.GetProp('_Name')
                new_name = f"{old_name}_rgrp-{i}"
                new_mol = Chem.MolFromSmiles(expansion)
                new_mol = add_props_to_new_mol(mol, new_mol)
                new_mol.SetProp('r_group', r_group)
                print(cmpd_id)
                new_mol.SetProp('cmpd_id', cmpd_id)
                new_mol.SetProp('smiles', expansion)
                new_mol.SetProp('smiles_before_rgrp_expansions', mol.GetProp('smiles'))
                new_mol.SetProp('_Name', new_name)
                new_mols.append(new_mol)

    print(len(new_mols), 'new mols found')
    end = time.time()
    w = Chem.SDWriter(os.path.join(args.output_dir, f'{args.name_tag}_r_group_expansions.sdf'))
    for new in new_mols:
        w.write(new)
    w.close()
    print('SDF file written')

    time_taken = round(end-start, 2)
    time_data = {'total_time': time_taken,
                 'n_cpus': args.n_cpus}
    dump_json(time_data, timings_file)


if __name__ == "__main__":
    main()
