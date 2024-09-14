import FragmentKnitwork.utils.fragmentConfig as config
import numpy as np
from FragmentKnitwork.utils.utils import get_distance
from rdkit import Chem
from rdkit.Chem import Mol, rdFMCS, rdShapeHelpers


def check_fragment_overlap(molA: Mol, molB: Mol, fragment_overlap=config.FRAGMENT_OVERLAP,
                           substructure_overlap=config.SUBSTRUCTURE_OVERLAP, mode='max'):
    """
    Check the % vol overlapping between fragments or substructures

    :param molA: rdkit mol
    :param molB: rdkit mol
    :param fragment_overlap: % volume overlap between fragments
    :param substructure_overlap: % volume overlap between substructures
    :param mode: whether to use the maximum overlap or the mean overlap as a threshold
    :return: whether molecule passes
    """
    if not molA or not molB:
        return False

    if mode not in ['max', 'mean']:
        raise ValueError("Mode must be max or mean")
    overlapA = 1 - rdShapeHelpers.ShapeProtrudeDist(molA, molB, allowReordering=False)
    overlapB = 1 - rdShapeHelpers.ShapeProtrudeDist(molB, molA, allowReordering=False)

    if mode == 'max':
        if max([overlapA, overlapB]) > fragment_overlap:
            return False
        else:
            return True

    if mode == 'mean':
        mean = (overlapA * overlapB) ** 0.5
        if mean > substructure_overlap:
            return False
        else:
            return True


def check_min_fragment_dist(molA: Mol, molB: Mol, distance: float = config.FRAGMENT_DISTANCE, return_distance: bool = False):
    """
    Check the min distance between fragments (by checking closest pair of atoms)

    :param molA:  rdkit mol
    :param molB: rdkit mol
    :param distance: maximum distance tolerated between fragments
    :return: whether molecule passes OR the min distance
    """
    if not molA or not molB:
        return False

    confA, confB = molA.GetConformer(), molB.GetConformer()
    distances = []

    for i in range(confA.GetNumAtoms()):
        posA = np.array(confA.GetAtomPosition(i))
        for j in range(confB.GetNumAtoms()):
            posB = np.array(confB.GetAtomPosition(j))
            distance = get_distance(posA, posB)
            distances.append(distance)

    if return_distance:
        return min(distances)

    if min(distances) > distance:
        return False
    else:
        return True


def get_num_carbons(smiles: str) -> int:
    """
    Count the number of carbons in a molecule

    :param smiles:
    :return:
    """
    c = Chem.MolFromSmarts("[#6]")
    return len(Chem.MolFromSmiles(smiles).GetSubstructMatches(c))


def run_r_group_expansion(smiles: str) -> dict:
    """
    Retrieves possible R groups for a molecule (for expansion purposes later) (a bit complicated but in the decomposition
    retrieves things that are removed before the constituent synthons are removed - would require 2 hops in transformation)

    :param smiles:
    :return: dict of synthon key and R groups that attach to that synthon
    """
    from FragmentKnitwork.utils.dbUtils import get_R_groups
    _synthons, _r_groups = get_R_groups(smiles)
    synthons = [syn for syn, r_group in zip(_synthons, _r_groups) if syn.count('[Xe]') == 1 and get_num_carbons(syn) >= 3]
    r_groups = [r_group for syn, r_group in zip(_synthons, _r_groups) if syn.count('[Xe]') == 1 and get_num_carbons(syn) >= 3]

    r_group_data = {synthon: [] for synthon in synthons}

    for synthon, r_group in zip(synthons, r_groups):
        # check the R groups are small so wouldn't have been considered in synthons
        if get_num_carbons(r_group) < 3 and '[Xe]' in r_group:
            if r_group not in r_group_data[synthon]:
                r_group_data[synthon].append(r_group)

    return r_group_data


def check_num_carbons(smiles: list, num_carbons: int = config.MIN_CARBONS) -> list:
    """
    Check the number of carbons in a substructure for merging

    :param smiles: list of substructure SMILES
    :param num_carbons: minimum number of carbons required
    :return: filtered list of smiles
    """
    c = Chem.MolFromSmarts("[#6]")
    filtered_smiles = [smi for smi in smiles if len(Chem.MolFromSmiles(smi).GetSubstructMatches(c)) >= num_carbons]
    return filtered_smiles


def check_carbon_ring(smiles: list, isSynthon: bool = False) -> list:
    """
    Check if molecule comprises a single carbon ring (for impure merging especially these queries tend to blow up)

    :param smiles: list of SMILES of molecules
    :param isSynthon: whether the molecule is a synthon or not (so has an attachment point)
    :return: list of filtered SMILES
    """
    filtered_smiles = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)

        num_rings = mol.GetRingInfo().NumRings()
        if num_rings != 1:
            filtered_smiles.append(smi)

        num_ring_atoms = sum([atom.IsInRing() for atom in mol.GetAtoms()])
        c = Chem.MolFromSmarts("[#6]")
        num_carbons = len(mol.GetSubstructMatches(c))
        num_atoms = mol.GetNumAtoms()

        if isSynthon:  # one atom will be a xenon (attachment point)
            if num_ring_atoms == (num_atoms - 1) and num_carbons == (num_atoms - 1):
                pass
            else:
                filtered_smiles.append(smi)

        else:
            if num_ring_atoms == num_atoms and num_carbons == num_atoms:
                pass
            else:
                filtered_smiles.append(smi)
    return filtered_smiles


def check_single_mol(smiles: list) -> list:
    """
    Check that the node is representing a single molecule (separate mols are denoted with full stop)

    :param smiles: list of SMILES of molecules
    :return: filtered list of SMILES
    """
    filtered_smiles = [smi for smi in smiles if '.' not in smi]
    return filtered_smiles


def get_mcs(larger_mol: Mol, sub_smi: Mol, return_mcs: bool = False):
    """
    Get the matching atoms between a larger molecule and a substructure (optionally return the MCS molecule)

    :param larger_mol: rdkit Mol for larger mol
    :param sub_smi:
    :param return_mcs: whether to return the mcs molecule
    :return: list of substructure matches (and optionally the mcs)
    """
    sub_mol = Chem.MolFromSmiles(sub_smi)
    mcs = Chem.MolFromSmarts(rdFMCS.FindMCS([larger_mol, sub_mol]).smartsString)
    matches = larger_mol.GetSubstructMatches(mcs)
    if return_mcs:
        return matches, mcs
    else:
        return matches


def atom_coords_from_pymol(mol_file: str) -> dict:
    """
    Get the coordinates and IDs of all atoms in a molfile

    :param mol_file:
    :return: dictionary with "coords" and "IDs"
    """
    from pymol import cmd
    coords = {"coords": [], "IDs": []}
    cmd.reinitialize()
    cmd.load(mol_file, "mol")
    cmd.iterate_state(1, "all", "coords.append([x,y,z])", space=coords)
    cmd.iterate_state(1, "all", "IDs.append(ID)", space=coords)
    return coords


def atom_IDs_to_molfile(mol_file, output_file, IDs):
    """
    Save a substructure of a mol file to a new mol files using the atom IDs

    :param mol_file:
    :param output_file:
    :param IDs:
    :return:
    """
    from pymol import cmd
    cmd.reinitialize()
    cmd.load(mol_file, "mol")
    ids_string = "id "
    for num in IDs[: len(IDs) - 1]:
        ids_string += str(num)
        ids_string += "+"
    ids_string += str(IDs[-1])
    cmd.select("substructure", ids_string)
    cmd.save(output_file, "substructure")
    return output_file


def round_coords(coords):
    return [round(c, 0) for c in coords]


def filt_by_bool(to_filter, bool_list):
    return [i for i, b in zip(to_filter, bool_list) if b]
