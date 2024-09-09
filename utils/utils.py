import json
import os

import FragmentKnitwork.utils.Config as config
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdmolfiles, rdShapeHelpers


def get_exps(data):
    all_exps = []
    for sub_pair in data:
        all_exps.extend(data[sub_pair]['expansions'])
    return all_exps

def load_json(fname):
    with open(fname, "r") as f:
        data = json.load(f)
    return data


def dump_json(data, fname):
    with open(fname, 'w') as f:
        json.dump(data, f)


def get_smiles(target, fragment, fragalysis_dir=config.FRAGALYSIS_DATA_DIR, fragalysis_version=config.FRAGALYSIS_VERSION):
    if fragalysis_version == 'v1':  # fpath format changed with new version of fragalysis download
        dir = os.path.join(fragalysis_dir, target, "aligned")
        fname_part = f"{target}-{fragment}"  # the first part of the filenames/folder names
        path = os.path.join(dir, fname_part)
        smiles_path = os.path.join(path, f"{fname_part}_smiles.txt")

    if fragalysis_version == 'v2':
        dir = os.path.join(fragalysis_dir, target, "aligned_files")
        path = os.path.join(dir, fragment)
        smiles_path = os.path.join(path, f"{fragment}_ligand.smi")

    try:
        with open(smiles_path) as smiles_file:  # open the file to get smiles
            smiles = smiles_file.read()
    except OSError as e:
        print("SMILES file cannot be found for that fragment.")
        print("Ensure files are saved correctly in Fragalysis format.")
        print(e)

    # smiles does not necessarily match what is in the network
    # get canonical smiles by converting to mol and converting back to smiles
    mol = Chem.MolFromSmiles(smiles)
    Chem.RemoveStereochemistry(mol)
    smiles = Chem.MolToSmiles(mol)

    return smiles


def get_mol(target, fragment, return_mol=False, fragalysis_dir=config.FRAGALYSIS_DATA_DIR,
            fragalysis_version=config.FRAGALYSIS_VERSION):
    if fragalysis_version == 'v1':
        dir = os.path.join(fragalysis_dir, target, "aligned")
        fname_part = f"{target}-{fragment}"  # the first part of the filenames/folder names
        path = os.path.join(dir, fname_part)
        mol_path = os.path.join(path, f"{fname_part}.mol")
    if fragalysis_version == 'v2':
        dir = os.path.join(fragalysis_dir, target, "aligned_files")
        path = os.path.join(dir, fragment)
        mol_path = os.path.join(path, f"{fragment}_ligand.mol")

    if return_mol:
        try:
            mol = rdmolfiles.MolFromMolFile(mol_path, sanitize=True)
            return mol
        except OSError as e:
            print("Mol file cannot be found for that fragment.")
            print("Ensure files are saved correctly in Fragalysis format.")
            print(e)
    else:
        return mol_path


def get_protein(target, fragment, return_mol=False, fragalysis_dir=config.FRAGALYSIS_DATA_DIR, protonated=False,
                desolv=True, fragalysis_version=config.FRAGALYSIS_VERSION):
    if fragalysis_version == 'v1':
        dir = os.path.join(fragalysis_dir, target, "aligned")
        fname_part = f"{target}-{fragment}"  # the first part of the filenames/folder names
    if fragalysis_version == 'v2':
        dir = os.path.join(fragalysis_dir, target, "aligned_files")
        fname_part = fragment  # the first part of the filenames/folder names
    if not protonated and desolv:
        protein_path = os.path.join(dir, fname_part, f"{fname_part}_apo-desolv.pdb")
    if protonated and desolv:
        protein_path = os.path.join(dir, fname_part, f"{fname_part}_apo-desolv-Hs.pdb")
    if not protonated and not desolv:
        protein_path = os.path.join(dir, fname_part, f"{fname_part}_apo.pdb")
    if protonated and not desolv:
        protein_path = os.path.join(dir, fname_part, f"{fname_part}_apo-Hs.pdb")

    if return_mol:
        try:
            mol = rdmolfiles.MolFromPDBFile(protein_path)
            return mol
        except OSError as e:
            print("Mol file cannot be found for that fragment.")
            print("Ensure files are saved correctly in Fragalysis format.")
            print(e)
    else:
        return protein_path


def get_bound_protein(
    target, fragment, return_mol=False, fragalysis_dir=config.FRAGALYSIS_DATA_DIR
):
    # bound-desolv files only used in fragalysis format v1
    from pymol import cmd
    dir = os.path.join(fragalysis_dir, target, "aligned")
    fname_part = f"{target}-{fragment}"  # the first part of the filenames/folder names
    protein_path = os.path.join(dir, fname_part, f"{fname_part}_bound-desolv.pdb")

    if os.path.exists(protein_path):
        if return_mol:
            mol = rdmolfiles.MolFromPDBFile(protein_path)
            return mol
        else:
            return protein_path

    else:
        print("Mol file cannot be found for that fragment. Will try creating file.")
        apo_file = get_protein(target, fragment, fragalysis_dir=fragalysis_dir)
        fragment_file = get_mol(target, fragment, fragalysis_dir=fragalysis_dir)
        if os.path.exists(apo_file) and os.path.exists(fragment_file):

            cmd.reinitialize()
            cmd.load(apo_file, "protein")
            cmd.load(fragment_file, "ligand")
            cmd.create('complex', 'ligand, protein')
            cmd.save(protein_path, "complex")
            mol = rdmolfiles.MolFromPDBFile(protein_path)

            if return_mol:
                mol = rdmolfiles.MolFromPDBFile(protein_path)
                return mol

            else:
                return protein_path

        else:
            print("Ensure files are saved correctly in Fragalysis format.")


def split_complex_file(holo_file, new_file):
    from pymol import cmd
    cmd.reinitialize()
    cmd.load(holo_file, 'complex')
    # cmd.select('lig', 'resn LIG')
    cmd.extract('hets', 'complex and HETATM')
    cmd.save(new_file, 'complex')
    return new_file

def get_distance(coord1, coord2):
    """
    :param coord1: 3D atom coordinates
    :param coord2: 3D atom coordinates
    :return: distance between the coordinates
    """
    sq = (coord1 - coord2) ** 2
    return np.sqrt(np.sum(sq))


def get_intersect(lst1, lst2):
    return list(set(lst1) & set(lst2))


def unnest_list(nested_list):
    return [i for sublist in nested_list for i in sublist]


def calculate_overlap(mol, ref_mol1):
    overlapA = 1 - rdShapeHelpers.ShapeProtrudeDist(ref_mol1, mol, allowReordering=False)
    return overlapA


def geom_mean(valA, valB):
    return (valA * valB) ** 0.5


def disable_rdlogger():
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')


def get_prop_as_list(mols, propName):
    return [mol.GetProp(propName) for mol in mols]


def order_by_lst(to_order, order_by, reverse=True):
    return [x for (y, x) in sorted(zip(order_by, to_order), reverse=reverse, key=lambda pair: pair[0])]