
import os
from argparse import ArgumentParser

from FragmentKnitwork.utils.quilterUtils import (add_props_from_dict,
                                                 get_ref_substructures)
from FragmentKnitwork.utils.utils import disable_rdlogger, load_json
from FragmentKnitwork.Quilter.scoreAlignedMols import process_molecule
from joblib import Parallel, delayed
from rdkit import Chem

disable_rdlogger()


def main():
    parser = ArgumentParser()
    parser.add_argument('--dir')
    parser.add_argument('--target')
    parser.add_argument('--n_cpus', type=int)
    parser.add_argument('--n_write', type=int, default=None)
    parser.add_argument('--output_dir')
    parser.add_argument('--substructure_dir')
    parser.add_argument('--diversity_selection', action='store_true')
    parser.add_argument('--type_merge', choices=['pure_merge', 'impure_merge', 'reverse_merge'])
    parser.add_argument('--butina_cluster_threshold', default=0.3)
    args = parser.parse_args()

    fragment_ints_file = os.path.join(args.substructure_dir, 'fragment_interactions.json')
    fragment_ints_data = load_json(fragment_ints_file)

    print('Reading mols')
    all_passing_mols, all_fA_mols, all_fB_mols = [], [], []

    pair_data = load_json(os.path.join(args.dir, "passing_data.json"))
    passing_mols, holo_files, json_files, mol_files, interaction_files, ref_subnode_mols, ref_synthon_mols = [], [], [], [], [], [], []

    # load data from saved json file
    names = pair_data['names']
    smiles = pair_data['smiles']
    cmpd_ids = pair_data['og_cmpd_id']
    ref_subnodes = pair_data['ref_subnode']
    ref_synthons = pair_data['ref_synthon']
    passes = pair_data['passes']
    if args.type_merge != 'pure_merge':
        similarities = pair_data['similarity']
        used_subnodes = pair_data['used_subnode']
        used_synthons = pair_data['used_synthon']
    else:
        similarities = ['NA'] * len(names)
        used_subnodes = ['NA'] * len(names)
        used_synthons = ['NA'] * len(names)
    if args.type_merge == 'reverse_merge':
        rev_synthons = pair_data['rev_synthon']
        intermed_nodes1 = pair_data['intermed_node1']
    else:
        rev_synthons = ['NA'] * len(names)
        intermed_nodes1 = ['NA'] * len(names)

    _fAs = pair_data['fragmentA']
    _fBs = pair_data['fragmentB']

    res_dirs = [os.path.join(args.dir, name.replace('_', '-')) for name in names]

    for res_dir, name, sim, smi, cmpd_id, used_subnode, used_synthon, ref_subnode, ref_synthon, \
        pas, rev_synthon, intermed_node, fA, fB in zip(res_dirs, names, similarities, smiles, \
        cmpd_ids, used_subnodes, used_synthons, ref_subnodes, ref_synthons, passes, rev_synthons, \
        intermed_nodes1, _fAs, _fBs):

        if pas:
            name_hyph = name.replace('_', '-')
            mol_file = os.path.join(res_dir, f"{name_hyph}.minimised.mol")
            mol = Chem.MolFromMolFile(mol_file)

            # add properties to the molecule
            if mol:
                mol.SetProp('_Name', name)
                mol.SetProp('similarity', str(sim))
                mol.SetProp('smiles', smi)
                mol.SetProp('cmpd_id', cmpd_id)
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
                holo_file = os.path.join(res_dir, f"{name_hyph}.holo_minimised.pdb")
                holo_files.append(holo_file)
                json_file = os.path.join(res_dir, f"{name_hyph}.minimised.json")
                json_files.append(json_file)
                interaction_file = os.path.join(res_dir, f"{name_hyph}_interactions.json")
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
            mol, fA=fA, fB=fB, substructsA=_ref_subnode_mols, substructsB=_ref_synthon_mols, target=args.target, \
            holo_file=holo_file, interactions_file=ints_file, fragment_ints_data=fragment_ints_data, mol_file=mol_file, \
            minimised_json_file=json_file) for mol, fA, fB, _ref_subnode_mols, _ref_synthon_mols, \
            holo_file, ints_file, mol_file, json_file in \
        zip(passing_mols, fAs, fBs, ref_subnode_mols, ref_synthon_mols, holo_files, interaction_files, mol_files, json_files)
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
