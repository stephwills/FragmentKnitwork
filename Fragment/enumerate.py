import itertools
import os

import FragmentKnitwork.utils.fragmentConfig as config
from FragmentKnitwork.Fragment.substructure import process_substructures_into_files
from FragmentKnitwork.utils.fragmentUtils import (check_fragment_overlap,
                                                  check_min_fragment_dist,
                                                  filt_by_bool)
from FragmentKnitwork.utils.utils import dump_json, load_json
from rdkit import Chem
from rdkit.Chem import Mol


def check_fragment_compatability(molA: Mol, molB: Mol, fragment_overlap=config.FRAGMENT_OVERLAP,
                                 distance=config.FRAGMENT_DISTANCE):
    """

    :param molA:
    :param molB:
    :param fragment_overlap:
    :param distance:
    :return:
    """
    # TODO: not currently used
    overlap_res = check_fragment_overlap(molA, molB, fragment_overlap, mode='max')
    if not overlap_res:
        return False

    dist_res = check_min_fragment_dist(molA, molB, distance)
    if not dist_res:
        return False

    return True


def generate_substructure_pairs(name_pair, pdbfile_pair, substructure_data, substructure_dir, target, prolif_check=config.PROLIF_CHECK,
                                substructure_overlap_check=config.SUBSTRUCTURE_OVERLAP_CHECK):
    """
    For a given pair of fragments, generate compatible pairs of substructures (depending on overlap)

    :param name_pair: tuple for the pair of fragments (named)
    :param pdbfile_pair: tuple with names of the corresponding pdb files
    :param substructure_data: dict containing the substructure data (fnames where they are saved)
    :param substructure_dir: where they are all saved
    :return: list of passing substructure pairs
    """
    fragmentA, fragmentB = name_pair[0], name_pair[1]
    protA, protB = pdbfile_pair[0], pdbfile_pair[1]
    subnodes = substructure_data[fragmentA]['subnode']
    synthons = substructure_data[fragmentB]['synthon']

    substructure_pairs = []

    num_substructure_pairs = 0
    num_good_overlap = 0
    num_good_ints = 0

    for subnode in subnodes:
        for synthon in synthons:
            num_substructure_pairs += 1

            # load the mols for the substructures
            subnode_files = subnodes[subnode]
            synthon_files = synthons[synthon]
            subnode_mols = [Chem.MolFromMolFile(os.path.join(substructure_dir, file)) for file in subnode_files]
            synthon_mols = [Chem.MolFromMolFile(os.path.join(substructure_dir, file)) for file in synthon_files]

            # check the degree of overlap between the substructures
            sub_pairs = list(itertools.product(subnode_mols, synthon_mols))
            # check if any of the possible substructure matches (may be multiple) have passing overlap
            if substructure_overlap_check:
                res = [check_fragment_overlap(pair[0], pair[1], mode='mean') for pair in sub_pairs]
                if sum(res) > 0:
                    num_good_overlap += 1
                    if prolif_check:
                        # check there is a combination of matches where both make interactions
                        from FragmentKnitwork.utils.quilterUtils import calc_prolif_interaction
                        res = []
                        for sub_pair in sub_pairs:
                            subnode_check = len(calc_prolif_interaction(sub_pair[0], protA, write_interactions_file=False, ligAsMol=True)) > 0
                            synthon_check = len(calc_prolif_interaction(sub_pair[1], protB, write_interactions_file=False, ligAsMol=True)) > 0
                            r = subnode_check and synthon_check
                            res.append(r)
                        if sum(res) > 0:
                            num_good_ints += 1
                            substructure_pairs.append([subnode, synthon])
                    else:
                        substructure_pairs.append([subnode, synthon])

            else:
                substructure_pairs.append([subnode, synthon])


    print('Number substructure pairs total:', num_substructure_pairs)
    print('Number substructure pairs after removing bad overlap:', num_good_overlap)
    print('Number substructure pairs after removing no ints:', num_good_ints)
    print('\n')
    return substructure_pairs


def enumeration(fragment_names: list,
                fragment_smiles: list,
                target: str,
                mols=None,
                mol_files=None,
                pdb_files=None,
                substructure_dir=None,
                overlap_check=config.CHECK_OVERLAP,
                distance_check=config.CHECK_DISTANCE,
                ignore_pairs=None):
    """
    Enumerate all pairs of substructures (and fragments) ready for querying the Fragment Network

    :param fragment_names: list of fragment names in form x1302_0A (number and chain)
    :param fragment_smiles: list of fragment smiles
    :param target: target (used for naming files)
    :param mols: list of fragment mols
    :param mol_files: list of mol files for the fragments
    :param substructure_dir: substructure dir to save all the enumerated substructures
    :param overlap_check: whether to check overlap
    :param distance_check: whether to check distance
    :return:
    """
    if not mols and not mol_files:
        raise ValueError('Provide mols or mol files for the fragments')

    if not mols:
        mols = [Chem.MolFromMolFile(file) for file in mol_files]

    # get pairs
    name_pairs = list(itertools.permutations(fragment_names, 2))
    smiles_pairs = list(itertools.permutations(fragment_smiles, 2))
    mol_pairs = list(itertools.permutations(mols, 2))
    pdbfile_pairs = list(itertools.permutations(pdb_files, 2))
    print(len(name_pairs), 'pairs in total')

    # if there are pairs to ignore, rule them out
    if ignore_pairs:
        filt_name_pairs = []
        filt_smiles_pairs = []
        filt_mol_pairs = []
        filt_pdbfile_pairs = []
        for name_pair, smiles_pair, mol_pair, pdb_pair in zip(name_pairs, smiles_pairs, mol_pairs, pdbfile_pairs):
            if list(name_pair) not in ignore_pairs:
                filt_name_pairs.append(name_pair)
                filt_smiles_pairs.append(smiles_pair)
                filt_mol_pairs.append(mol_pair)
                filt_pdbfile_pairs.append(pdb_pair)
        print(len(filt_name_pairs), 'remaining after removing', len(ignore_pairs), 'pairs to ignore')
        name_pairs = filt_name_pairs
        smiles_pairs = filt_smiles_pairs
        mol_pairs = filt_mol_pairs
        pdbfile_pairs = filt_pdbfile_pairs

    # TODO: or use single function above
    # check overlap and distance between the fragments in the pairs
    if overlap_check:
        res = [check_fragment_overlap(pair[0], pair[1]) for pair in mol_pairs]
        mol_pairs = filt_by_bool(mol_pairs, res)
        name_pairs = filt_by_bool(name_pairs, res)
        smiles_pairs = filt_by_bool(smiles_pairs, res)
        pdbfile_pairs = filt_by_bool(pdbfile_pairs, res)
        print(len(name_pairs), 'pairs after removing overlapping mols')

    if distance_check:
        res = [check_min_fragment_dist(pair[0], pair[1]) for pair in mol_pairs]
        mol_pairs = filt_by_bool(mol_pairs, res)
        name_pairs = filt_by_bool(name_pairs, res)
        smiles_pairs = filt_by_bool(smiles_pairs, res)
        pdbfile_pairs = filt_by_bool(pdbfile_pairs, res)
        print(len(name_pairs), 'pairs after removing far apart mols')

    # generate substructures if not generated already (checks if directory has any files)
    substructure_data_fname = os.path.join(substructure_dir, 'substructure_files.json')
    if os.path.exists(substructure_data_fname):
        print('Substructure files already in substructure dir -- not regenerating, please check')

    else:
        print('Generating substructure files')
        process_substructures_into_files(fragment_names,
                                         fragment_smiles,
                                         mol_files,
                                         substructure_dir)
        print('Substructure files generated')

    substructure_data = load_json(substructure_data_fname)
    substructure_pair_data = {}

    for name_pair, smiles_pair, pdbfile_pair in zip(name_pairs, smiles_pairs, pdbfile_pairs):
        substructure_pairs = generate_substructure_pairs(name_pair,
                                                         pdbfile_pair,
                                                         substructure_data,
                                                         substructure_dir,
                                                         target)
        if len(substructure_pairs) > 0:
            pair_name = f"{name_pair[0]}-{name_pair[1]}"
            substructure_pair_data[pair_name] = substructure_pairs

    print(len(substructure_pair_data), 'fragment pairs with expansions after enumeration')
    c = 0
    for pair in substructure_pair_data:
        c += len(substructure_pair_data[pair])
    print(c, 'substructure pairs for querying after enumeration')
    substructure_pair_fname = os.path.join(substructure_dir, 'substructure_pairs.json')
    dump_json(substructure_pair_data, substructure_pair_fname)
    return substructure_pair_fname
