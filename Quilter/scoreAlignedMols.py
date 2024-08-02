
import os
from argparse import ArgumentParser

from FragmentKnitwork.Quilter.checkPlacedMol import process_molecule
from FragmentKnitwork.utils.quilterUtils import (add_props_from_dict,
                                                 get_ref_substructures)
from FragmentKnitwork.utils.utils import disable_rdlogger, load_json
from joblib import Parallel, delayed
from rdkit import Chem
from tqdm import tqdm

disable_rdlogger()


def main():
    """
    Currently doing this because ProLIF not working in the same environment used for conformer generation

    :return:
    """
    parser = ArgumentParser()
    parser.add_argument('--dir')
    parser.add_argument('--target')
    parser.add_argument('--n_cpus', type=int)
    parser.add_argument('--n_write', type=int, default=None)
    parser.add_argument('--output_dir')
    parser.add_argument('--substructure_dir')
    parser.add_argument('--type_merge', choices=['pure_merge', 'impure_merge', 'reverse_merge'])
    args = parser.parse_args()

    fragment_ints_data = load_json(os.path.join(args.substructure_dir, 'fragment_interactions.json'))

    print('Reading mols')
    all_passing_mols, all_fA_mols, all_fB_mols = [], [], []

    # read in all the pairs from the ouptut directory
    pair_names = [i for i in os.listdir(args.dir) if i[0] == 'x']
    for pair_name in tqdm(pair_names, position=0, leave=True, total=len(pair_names)):
        print('Running', pair_name)
        pair_dir = os.path.join(args.dir, pair_name)

        fA, fB = pair_name.split('-')[0], pair_name.split('-')[1]
        pair_data = load_json(os.path.join(pair_dir, f"passing_data.json"))
        passing_mols, holo_files, json_files, mol_files, interaction_files, ref_subnode_mols, ref_synthon_mols = [], [], [], [], [], [], []

        # load data from saved json file
        names = pair_data['names']
        smiles = pair_data['smiles']
        cmpd_ids = pair_data['cmpd_ids']
        ref_subnodes = pair_data['ref_subnodes']
        ref_synthons = pair_data['ref_synthons']
        passes = pair_data['passes']
        if args.type_merge != 'pure_merge':
            if 'similarities' in pair_data.keys():
                similarities = pair_data['similarities']
            else:
                similarities = ['NA'] * len(names)
            used_subnodes = pair_data['used_subnodes']
            used_synthons = pair_data['used_synthons']
        else:
            similarities = ['NA'] * len(names)
            used_subnodes = ['NA'] * len(names)
            used_synthons = ['NA'] * len(names)
        if args.type_merge == 'reverse_merge':
            rev_synthons = pair_data['rev_synthons']
            intermed_nodes1 = pair_data['intermed_nodes1']
        else:
            rev_synthons = ['NA'] * len(names)
            intermed_nodes1 = ['NA'] * len(names)

        for name, sim, smi, cmpd_id, used_subnode, used_synthon, \
            ref_subnode, ref_synthon, pas, rev_synthon, intermed_node \
                in zip(names, similarities, smiles, cmpd_ids, used_subnodes, used_synthons, \
                       ref_subnodes, ref_synthons, passes, rev_synthons, intermed_nodes1):
            if pas:
                name_hyph = name.replace('_', '-')
                mol_file = os.path.join(pair_dir, name_hyph, f"{name_hyph}.minimised.mol")
                mol = Chem.MolFromMolFile(os.path.join(pair_dir, name_hyph, f"{name_hyph}.minimised.mol"))
                if mol:

                    # add properties to the molecule
                    mol.SetProp('_Name', name)
                    mol.SetProp('similarity', str(sim))
                    mol.SetProp('smiles', smi)
                    mol.SetProp('cmpd_id', ','.join(cmpd_id))
                    mol.SetProp('used_subnode', used_subnode)
                    mol.SetProp('used_synthon', used_synthon)
                    mol.SetProp('ref_subnode', ref_subnode)
                    mol.SetProp('ref_synthon', ref_synthon)
                    mol.SetProp('fragmentA', fA)
                    mol.SetProp('fragmentB', fB)
                    if args.type_merge == 'reverse_merge':
                        mol.SetProp('rev_synthon', rev_synthon)
                        mol.SetProp('intermed_node1', intermed_node)
                    passing_mols.append(mol)
                    name_hyph = name.replace('_', '-')

                    # get file names for processing
                    holo_file = os.path.join(pair_dir, name_hyph, f"{name_hyph}.holo_minimised.pdb")
                    holo_files.append(holo_file)
                    json_file = os.path.join(pair_dir, name_hyph, f"{name_hyph}.minimised.json")
                    json_files.append(json_file)
                    interaction_file = os.path.join(pair_dir, name_hyph, f"{name_hyph}_interactions.json")
                    interaction_files.append(interaction_file)
                    mol_files.append(mol_file)

                    # load substructures used for placement
                    if args.type_merge == 'impure_merge' or args.type_merge == 'pure_merge':
                        _ref_subnode_mols = get_ref_substructures(ref_subnode, fA, substructure_dir=args.substructure_dir, isSynthon=False)
                        _ref_synthon_mols = get_ref_substructures(ref_synthon, fB, substructure_dir=args.substructure_dir, isSynthon=True)
                    else:
                        _ref_subnode_mols = get_ref_substructures(ref_subnode, fA, substructure_dir=args.substructure_dir, isSynthon=True)
                        _ref_synthon_mols = get_ref_substructures(ref_synthon, fB, substructure_dir=args.substructure_dir, isSynthon=True)
                    ref_subnode_mols.append(_ref_subnode_mols)
                    ref_synthon_mols.append(_ref_synthon_mols)

        print('Read', len(passing_mols), 'mols')
        fAs = [mol.GetProp('fragmentA') for mol in passing_mols]
        fBs = [mol.GetProp('fragmentB') for mol in passing_mols]

        print('Calculating checks')
        results = Parallel(n_jobs=args.n_cpus, backend='multiprocessing')(
            delayed(process_molecule)(
                mol, fA=fA, fB=fB, substructsA=_ref_subnode_mols, substructsB=_ref_synthon_mols, \
                target=args.target, holo_file=holo_file, interactions_file=ints_file, \
                fragment_ints_data=fragment_ints_data, mol_file=mol_file, minimised_json_file=json_file
            ) for mol, fA, fB, _ref_subnode_mols, _ref_synthon_mols, holo_file, ints_file, \
                  mol_file, json_file in zip(passing_mols, fAs, fBs, ref_subnode_mols, ref_synthon_mols, \
                                             holo_files, interaction_files, mol_files, json_files)
        )
        print('Calculated properties')

        # get only the passing mols
        passing_mols = [mol for mol, res in zip(passing_mols, results) if res[0]]
        results = [res for res in results if res[0]]
        data_dicts = [r[1] for r in results]
        fA_mols = [r[2] for r in results]
        fB_mols = [r[3] for r in results]

        for mol, data_dict, fA_mol, fB_mol in zip(passing_mols, data_dicts, fA_mols, fB_mols):
            new_mol = add_props_from_dict(mol, data_dict)
            all_passing_mols.append(new_mol)
            all_fA_mols.append(fA_mol)
            all_fB_mols.append(fB_mol)

    w = Chem.SDWriter(os.path.join(args.output_dir, 'ordered_mols.sdf'))
    wA = Chem.SDWriter(os.path.join(args.output_dir, 'ordered_fA_mols.sdf'))
    wB = Chem.SDWriter(os.path.join(args.output_dir, 'ordered_fB_mols.sdf'))

    for mol, fA_mol, fB_mol in zip(all_passing_mols, all_fA_mols, all_fB_mols):
        w.write(mol)
        wA.write(fA_mol)
        wB.write(fB_mol)

    print(len(all_passing_mols), 'mols')
    print('Files written')


if __name__ == "__main__":
    main()