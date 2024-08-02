
from FragmentKnitwork.Quilter.checkInteractions import check_interactions
from FragmentKnitwork.Quilter.checkOverlap import check_overlap
from FragmentKnitwork.utils.quilterUtils import check_energy, mean_SuCOS
from FragmentKnitwork.utils.utils import disable_rdlogger, load_json
from rdkit.Chem import Crippen, rdMolDescriptors

disable_rdlogger()


def calculate_descriptors(mol):
    """

    :param mol:
    :return:
    """
    mol_data = {}
    HA = mol.GetNumHeavyAtoms()
    mol_data['HA'] = HA
    mol_data['HBD'] = rdMolDescriptors.CalcNumHBD(mol)
    mol_data['HBA'] = rdMolDescriptors.CalcNumHBA(mol)
    mol_data['MW'] = rdMolDescriptors.CalcExactMolWt(mol)
    mol_data['TPSA'] = rdMolDescriptors.CalcTPSA(mol)
    mol_data['LogP'] = Crippen.MolLogP(mol)
    mol_data['RBs'] = rdMolDescriptors.CalcNumRotatableBonds(mol)
    return mol_data


def add_data_from_fragmenstein_json(mol_data, minimised_json_file):
    """

    :param mol_data:
    :param minimised_json_file:
    :return:
    """
    data = load_json(minimised_json_file)
    G_bound = data["Energy"]["bound"]["total_score"]  # energy bound
    G_unbound = data["Energy"]["unbound"]["total_score"]  # unbound
    deltaG = G_bound - G_unbound  # calculate energy difference
    comRMSD = data["mRMSD"]  # RMSD between two fragments and merge
    mol_data['G_bound'] = G_bound
    mol_data['G_unbound'] = G_unbound
    mol_data['deltaG'] = deltaG
    mol_data['comRMSD'] = comRMSD
    return mol_data


def process_molecule(mol, fA, fB, substructsA, substructsB, target, holo_file, interactions_file, fragment_ints_data,
                     mol_file, minimised_json_file):
    """

    :param mol:
    :param fA:
    :param fB:
    :param substructsA:
    :param substructsB:
    :param target:
    :param holo_file:
    :param interactions_file:
    :param fragment_ints_data:
    :param mol_file:
    :param minimised_json_file:
    :return:
    """
    mol_data = calculate_descriptors(mol)
    mol_data = add_data_from_fragmenstein_json(mol_data, minimised_json_file)

    overlap_check = check_overlap(mol, substructsA, substructsB)
    mol_data['substructure_overlap_check'] = overlap_check

    energy_check, energy_ratio = check_energy(mol, returnRatio=True)
    mol_data['energy_ratio_check'] = energy_check
    mol_data['const_energy/unconst_energy'] = energy_ratio

    sucos, fA_mol, fB_mol = mean_SuCOS(mol, fA=fA, fB=fB, target=target)
    sucos_check = sucos >= 0.55
    mol_data['sucos'] = sucos
    mol_data['sucos_check'] = sucos_check

    interaction_check, interaction_data = check_interactions(holo_file, mol_file, fA, fB, fragment_ints_data, interactions_file)
    mol_data.update(interaction_data)

    if energy_check and overlap_check:
        pass_checks = True
    else:
        pass_checks = False

    return pass_checks, mol_data, fA_mol, fB_mol
