import json
import os
from typing import List

from FragmentKnitwork.utils import fragmentConfig as config
from FragmentKnitwork.utils.dbUtils import get_subnodes, get_synthons
from FragmentKnitwork.utils.fragmentUtils import (atom_coords_from_pymol,
                                                  atom_IDs_to_molfile,
                                                  check_carbon_ring,
                                                  check_num_carbons,
                                                  check_single_mol, get_mcs,
                                                  round_coords)
from FragmentKnitwork.utils.utils import dump_json
from rdkit import Chem
from rdkit.Chem import Mol


def attachment_idxs(mol: Mol, substruct_matches: List) -> List:
    """
    Get the index of atoms attached to the substructure (could be multiple possible matches and multiple atoms)

    :param mol:
    :param substruct_matches:
    :return: list of list of attachment atom indices
    """
    # TODO: not currently used
    attachment_atoms = []

    for matches in substruct_matches:
        match_attachment_atoms = []
        for atom_idx in matches:
            for neigh in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                if neigh.GetIdx() not in matches:
                    match_attachment_atoms.append(atom_idx)
        attachment_atoms.append(match_attachment_atoms)

    return attachment_atoms


def get_substructure_atom_coords(mol_file: str, substructure: str):
    """
    Get the atom coordinates for the substructure from the mol file

    :param mol_file: mol file for the fragment
    :param substructure: SMILES of the substructure
    :return:
    """
    mol = Chem.MolFromMolFile(mol_file)
    substructure_matches, mcs = get_mcs(mol, substructure, return_mcs=True)

    if not substructure_matches or len(substructure_matches) == 0:
        print("Could not find MCS between the mol and the substructure")
        return None, None

    # TODO: get attachment points

    # for each possible substructure match, get the atom coordinates using RDKit
    substructure_atom_coords = []
    for matching_atoms in substructure_matches:
        coords = []
        for atom in matching_atoms:
            coord = list(mol.GetConformer().GetAtomPosition(atom))
            coords.append(coord)
        substructure_atom_coords.append(coords)

    return substructure_atom_coords


def save_substructure_to_file(substructure: str, mol_file: str, pymol_coords, type, fragment, dir, to_round=True):
    """
    Save a substructure to a new mol file by checking for correlation between RDKit and PyMOL coordinates
    (there are probably better ways to do this)

    :param substructure: smiles of the substructure
    :param mol_file: fragment mol file substructure coming from
    :param pymol_coords: list of coordinates for all atoms in the molecule
    :param type: whether a synthon or subnode
    :param fragment: fragment name
    :param dir: where to save the substructure file
    :param to_round:
    :return: list of output files (might be multiple matches)
    """
    substructure_atom_coords = get_substructure_atom_coords(mol_file,
                                                            substructure)

    output_files = []
    if to_round:
        rounded = [round_coords(coords) for coords in pymol_coords["coords"]]
        pymol_coords["coords"] = rounded

    if substructure_atom_coords: # [[coord, coord, coord], [coord, coord, coord]]
        for i, atom_coords in enumerate(substructure_atom_coords):
            substructure_IDs = []
            # print('atom coords', atom_coords)
            atom_coords = [round_coords(coords) for coords in atom_coords]
            for coords in atom_coords:
                if coords in pymol_coords["coords"]:  # get the pymol ID associated with the coordinate
                    idx = pymol_coords["coords"].index(coords)
                    ID = pymol_coords["IDs"][idx]
                    substructure_IDs.append(ID)

            fname = f"{type}_{fragment}_{substructure}_{i}.mol"
            fname = os.path.join(dir, fname)

            if len(substructure_IDs) > 0:  # use the atom IDs of the substructure to create new mol file with pymol
                output_file = atom_IDs_to_molfile(mol_file, fname, substructure_IDs)
                output_files.append(os.path.basename(output_file))

    return output_files


def equivalent_synthon_subnodes(fragment_names, fragment_smiles, returnMatches=True, saveJson=None):
    """
    Record the equivalent synthon for a given subnode (when we want to replace both substructures in a merge)

    :param fragment_names:
    :param fragment_smiles:
    :param returnMatches:
    :param saveJson:
    :return:
    """
    all_matches = {}

    for fragment_name, fragment_smiles in zip(fragment_names, fragment_smiles):
        synthon_smiles = get_synthons(fragment_smiles, terminal_synthons=False)
        subnode_smiles = get_subnodes(fragment_smiles, terminal_subnodes=False)

        synthons = [Chem.MolFromSmiles(syn) for syn in synthon_smiles]
        subnodes = [Chem.MolFromSmiles(sub) for sub in subnode_smiles]

        smiles_matches = {}
        for i, subnode in enumerate(subnodes):
            matches = []
            for j, synthon in enumerate(synthons):
                num_atoms = synthon.GetNumAtoms()
                num_match_atoms = synthon.GetSubstructMatch(subnode)
                if len(num_match_atoms) == num_atoms - 1:
                    matches.append(synthon_smiles[j])
            if len(matches) > 0:
                smiles_matches[subnode_smiles[i]] = ','.join(matches)
        all_matches[fragment_name] = smiles_matches

    if saveJson:
        dump_json(all_matches, saveJson)

    if returnMatches:
        return all_matches


def process_substructures_into_files(fragment_names: list, fragment_smiles: list, mol_files: list, output_dir: str,
                                     carbons_check: bool = config.CARBONS_CHECK, single_mol_check: bool = config.SINGLE_MOL_CHECK,
                                     carbon_ring_check: bool = config.SINGLE_CARBON_RING_CHECK,
                                     custom_subnodes=None, custom_synthons=None,
                                     terminal_subnodes: bool = config.TERMINAL_SUBNODES, terminal_synthons: bool = config.TERMINAL_SYNTHONS):
    """
    Given the fragments, process them into substructures using the Fragment Network and save the substructures
    with coordinates to mol files (may be multiple possible matches for each)

    :param fragment_names: list of fragment names
    :param fragment_smiles: list of fragment smiles
    :param mol_files: list of mol files associated with the fragments
    :param output_dir: output dir to save the substructures as files
    :param carbons_check: whether to check the number of carbons in the substructure for it to pass
    :param custom_subnodes: optional (provide custom subnodes)
    :param custom_synthons: optional (provide custom synthons)
    :return:
    """
    substructure_files = {}

    # enumerate substructures using the Fragment Network
    for fragment, smi, mol_file in zip(fragment_names, fragment_smiles, mol_files):
        subnodes = get_subnodes(smi, terminal_subnodes=terminal_subnodes)
        synthons = get_synthons(smi, terminal_synthons=terminal_synthons)

        if carbons_check:
            subnodes = check_num_carbons(subnodes)
            synthons = check_num_carbons(synthons)
        if single_mol_check:
            subnodes = check_single_mol(subnodes)
            synthons = check_single_mol(synthons)
        if carbon_ring_check:
            subnodes = check_carbon_ring(subnodes)
            synthons = check_carbon_ring(synthons, True)

        pymol_coords = atom_coords_from_pymol(mol_file)  # get the atomic coordinates using pymol API

        if custom_synthons:
            synthons = custom_synthons
        if custom_subnodes:
            subnodes = custom_subnodes

        # save possible substructures as individual files and store fnames in dict
        subnode_files = {}
        synthon_files = {}
        # create mol files for the matching substructures in a designated directory
        for substructure in subnodes:
            output_files = save_substructure_to_file(substructure, mol_file, pymol_coords,
                                                     'subnode', fragment, output_dir)
            subnode_files[substructure] = output_files
        for substructure in synthons:
            output_files = save_substructure_to_file(substructure, mol_file, pymol_coords,
                                                     'synthon', fragment, output_dir)
            synthon_files[substructure] = output_files

        substructure_files[fragment] = {'subnode': subnode_files,
                                        'synthon': synthon_files}

    data_fname = os.path.join(output_dir, 'substructure_files.json')
    dump_json(substructure_files, data_fname)
