import itertools
import os
import shutil

import numpy as np
from FragmentKnitwork.utils import quilterConfig as config
from FragmentKnitwork.utils.utils import (dump_json, geom_mean, get_mol,
                                          load_json)
from rdkit import Chem
from rdkit.Chem import AllChem, ChemicalFeatures, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps


fdef = AllChem.BuildFeatureFactory(config.ALIGNMENT_FDEF)
fmParams = {}
for k in fdef.GetFeatureFamilies():
    fparams = FeatMaps.FeatMapParams()
    fmParams[k] = fparams


def create_alignment_dirs(fragmentA, fragmentB, working_dir, output_dir):
    """
    Create individual pair dirs for fragment pair results

    :param fragmentA:
    :param fragmentB:
    :param working_dir:
    :param output_dir:
    :return:
    """
    pair = fragmentA + '-' + fragmentB
    pair_working_dir = os.path.join(working_dir, pair)
    pair_output_dir = os.path.join(output_dir, pair)
    if not os.path.exists(pair_working_dir):
        os.mkdir(pair_working_dir)
    if not os.path.exists(pair_output_dir):
        os.mkdir(pair_output_dir)
    return pair, pair_working_dir, pair_output_dir


def get_FeatureMapScore(small_m, large_m, score_mode=FeatMaps.FeatMapScoreMode.Best):
    """
    # Code from https://github.com/MarcMoesser/SuCOS

    :param small_m:
    :param large_m:
    :param score_mode:
    :param fdef:
    :return:
    """

    featLists = []
    for m in [small_m, large_m]:
        rawFeats = fdef.GetFeaturesForMol(m)
        featLists.append([f for f in rawFeats])

    fms = [FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params=fmParams) for x in featLists]
    fms[0].scoreMode = score_mode
    fm_score = fms[0].ScoreFeats(featLists[1]) / min(fms[0].GetNumFeatures(), len(featLists[1]))
    return fm_score


def calc_SuCOS(reflig, prb_mol, score_mode=FeatMaps.FeatMapScoreMode.Best, return_all=False):
    """
    # Code from https://github.com/MarcMoesser/SuCOS

    :param reflig:
    :param prb_mol:
    :param score_mode:
    :param return_all:
    :return:
    """
    fm_score = get_FeatureMapScore(reflig, prb_mol, score_mode)
    protrude_dist = rdShapeHelpers.ShapeProtrudeDist(reflig, prb_mol, allowReordering=False)
    SuCOS_score = 0.5 * fm_score + 0.5 * (1 - protrude_dist)

    if return_all:
        return SuCOS_score, fm_score, (1 - protrude_dist)
    else:
        return SuCOS_score


def mean_SuCOS(mol, fA_mol=None, fB_mol=None, fA=None, fB=None, target=None, asMol=False, returnFragmentMols=True):
    """
    Calculate the mean SuCOS score (shape and colour overlap) using RDKit against two reference molecules.

    :param mol:
    :param fA_mol:
    :param fB_mol:
    :param fA:
    :param fB:
    :param target:
    :param asMol:
    :param returnFragmentMols:
    :return:
    """
    if not asMol:
        fA_mol = get_mol(target, fA, True)
        fB_mol = get_mol(target, fB, True)
    sucosA = calc_SuCOS(fA_mol, mol)
    sucosB = calc_SuCOS(fB_mol, mol)
    mean = (sucosA * sucosB) ** 0.5
    if not returnFragmentMols:
        return mean
    else:
        return mean, fA_mol, fB_mol


def get_scoring_function(score=config.ALIGNMENT_SCORE):
    """
    Redundant at this time as we only have one scoring function implemented.

    :param score:
    :return:
    """
    if score == 'SuCOS':
        return calc_SuCOS


def eval_conformer(conformer, ref_mol1, ref_mol2, score=config.ALIGNMENT_SCORE):
    """
    Redundant at this time as we only have one scoring function implemented.

    :param conformer:
    :param ref_mol1:
    :param ref_mol2:
    :param score:
    :return:
    """
    scoring_function = get_scoring_function(score)
    return geom_mean(scoring_function(ref_mol1, conformer), scoring_function(ref_mol2, conformer))


def load_feat_factory(fdef=config.ALIGNMENT_FDEF):
    """
    Load feature factory using the alignment fdef

    :param fdef:
    :return:
    """
    featFactory = ChemicalFeatures.BuildFeatureFactory(fdef)
    return featFactory


def get_all_length_combinations(lst):
    """

    :param lst:
    :return:
    """
    combinations = []
    for i in range(len(lst)):
        combinations.extend([list(l) for l in (itertools.combinations(lst, i + 1))])
    return combinations


def all_combinations_pharmacophores(pharmacophores_A, pharmacophores_B):
    """
    Get all possible combinations of pharmacophore features for creating an embedding

    :param pharmacophores_A:
    :param pharmacophores_B:
    :return:
    """
    combinations_A = get_all_length_combinations(pharmacophores_A)
    combinations_B = get_all_length_combinations(pharmacophores_B)

    unique_combinations = [[i, j] for j in combinations_B for i in combinations_A]
    return unique_combinations


def get_best(items, scores, scoring_mode):
    """
    Get the best score depending on scoring mode

    :param items:
    :param scores:
    :param scoring_mode:
    :return:
    """
    if scoring_mode == 'min':
        best_score = min(scores)
    elif scoring_mode == 'max':
        best_score = max(scores)
    idx = scores.index(best_score)
    return best_score, idx, items[idx]


def loop_list(lst):
    """

    :param lst:
    :return:
    """
    mappings = []
    for i in range(len(lst)):
        mappings.append(lst[i:len(lst)] + lst[0:i])
    return mappings


def coords_from_ids(atom_ids, mol):
    """
    Get atom coords associated with specific atom IDs
    :param atom_ids:
    :param mol:
    :return:
    """
    return [np.array(mol.GetConformer().GetAtomPosition(idx)) for idx in atom_ids]


def get_ref_substructures(smiles, fragment, isSynthon=False, asMols=True, substructure_dir=config.SUBSTRUCTURE_DIR, info_file=config.INFO_FILE_EXISTS):
    """
    Retrieve the possible reference substructures extracted from the original fragment mols using the file structure

    :param smiles: smiles of the substrcuture
    :param fragment: name of the fragment
    :param isSynthon: whether a synthon or subnode (synthon has attachment point included)
    :param asMols: return as mols rather than fnames
    :param substructure_dir: where all the substructures have been saved
    :return: either list of mols or fnames
    """
    # load more quickly using a pre-generated json containing dict with files associated with each substructure SMILES
    if info_file:
        data = load_json(os.path.join(substructure_dir, 'substructure_files.json'))
        if isSynthon:
            fnames = data[fragment]['synthon'][smiles]
        else:
            fnames = data[fragment]['subnode'][smiles]
    # otherwise have to check file names
    else:
        all_files = os.listdir(substructure_dir)
        if isSynthon:
            fnames = [file for file in all_files if 'synthon' in file and smiles in file and fragment in file]
        else:
            fnames = [file for file in all_files if 'subnode' in file and smiles in file and fragment in file]
    fnames = [os.path.join(substructure_dir, file) for file in fnames]

    if asMols:
        # some of the molecules don't load and need rectifying
        mols = []
        for file in fnames:
            mol = Chem.MolFromMolFile(file)
            if mol:
                mols.append(mol)
            else:
                # if mol doesn't load, use molecular rectifier
                from molecular_rectifier import Rectifier
                mol = Chem.MolFromMolFile(file, sanitize=False)
                recto = Rectifier(mol)
                recto.fix()
                fixed_mol = recto.mol
                fixed_mol = Chem.RemoveHs(fixed_mol)
                mols.append(fixed_mol)

        return mols
    else:
        return fnames


def move_files(merge_dirs, output_dir=config.OUTPUT_DIR, minimal_files=config.MIN_FILES):
    """
    Move files from working dir to output dir (helps prevent from keeping all the Fragmenstein output files that aren't
    needed)
    :param merge_dirs:
    :param output_dir:
    :param minimal_files:
    :return:
    """
    if not minimal_files:
        for merge_dir in merge_dirs:
            shutil.move(merge_dir, os.path.join(output_dir, os.path.basename(merge_dir)))

    else:
        for merge_dir in merge_dirs:
            bname_dir = os.path.basename(merge_dir)
            new_merge_dir = os.path.join(output_dir, bname_dir)
            if not os.path.exists(new_merge_dir):
                os.mkdir(new_merge_dir)

            json_file = f"{bname_dir}.minimised.json"
            mol_file = f"{bname_dir}.minimised.mol"
            prot_file = f"{bname_dir}.holo_minimised.pdb"

            bname_files = [json_file, mol_file, prot_file]
            for file in bname_files:
                shutil.move(os.path.join(merge_dir, file), os.path.join(new_merge_dir, file))
            shutil.rmtree(merge_dir)


def add_props_to_new_mol(original_mol, new_mol, property_dict=None, name=None):
    """
    Function to add properties from a dictionary or an existing molecule to a new molecule (because properties not
    always preserved so have to be re-added)

    :param original_mol:
    :param new_mol:
    :param property_dict: dictionary to add properties from (instead of an old molecule)
    :param name: name of molecule in property dict if using
    :return: new mol with properties assigned
    """
    if property_dict:
        for prop in property_dict[name]:
            new_mol.SetProp(prop, property_dict[name][prop])
    else:
        prop_names = list(original_mol.GetPropNames())
        for prop_name in prop_names:
            prop = original_mol.GetProp(prop_name)
            new_mol.SetProp(prop_name, prop)
    return new_mol


def add_props_from_dict(mol, property_dict):
    """
    Add properties to molecule directly from dictionary

    :param mol:
    :param property_dict:
    :return:
    """
    for prop_name in property_dict:
        mol.SetProp(prop_name, str(property_dict[prop_name]))
    return mol


def check_for_overlap(mol, substructA, substructB):
    """
    Check that there is overlap between mol and the two substructures used for its placement (there should be)

    :param mol:
    :param substructA:
    :param substructB:
    :return:
    """
    overlapA = 1 - rdShapeHelpers.ShapeProtrudeDist(mol, substructA)
    overlapB = 1 - rdShapeHelpers.ShapeProtrudeDist(mol, substructB)
    if overlapA == 0 or overlapB == 0:
        return False
    else:
        return True


def calc_energy(mol):
    """
    Calculate energy of molecule using UFF

    :param mol:
    :return:
    """
    mol_energy = AllChem.UFFGetMoleculeForceField(mol).CalcEnergy()
    return mol_energy


def calc_unconstrained_energy(og_mol, n_conf=config.ENERGY_CHECK_N_CONF):
    """
    Calculate average energy of molecule based on several conformers

    :param og_mol:
    :param n_conf:
    :return:
    """
    unconstrained_energies = []
    for i in range(n_conf):  # generate conformations and calculate energy
        mol = Chem.Mol(og_mol)
        AllChem.EmbedMolecule(mol, randomSeed=i)
        AllChem.UFFOptimizeMolecule(mol)
        e = calc_energy(mol)
        unconstrained_energies.append(e)

    # calculate the average of all the energies
    avg = sum(unconstrained_energies) / len(unconstrained_energies)
    return avg


def check_energy(mol, n_conf: int = config.ENERGY_CHECK_N_CONF, energy_threshold: float = config.ENERGY_CHECK_THRESHOLD, returnRatio=False):
    """
    Run the energy check (checking energy of the aligned mol against the energies of several unconstrained
    conformations)

    :param mol: molecule to check
    :param n_conf: number of unconstrained conformations to generate
    :param energy_threshold:
    :return: whether to return the ratio
    """
    try:
        const_energy = calc_energy(mol)
        unconst_energy = calc_unconstrained_energy(mol, n_conf)
        ratio = const_energy / unconst_energy
        if ratio >= energy_threshold:
            result = False
        else:
            result = True
        if not returnRatio:
            return result
        else:
            return result, ratio
    except Exception as e:
        print(e)
        if not returnRatio:
            return False
        else:
            return False, None


def make_complex_file(new_fname, fragment_file: str, apo_file: str, dir: str) -> str:
    """
    Create pdb file containing complex of fragment and the apo minimized protein for PLIP

    :param new_fname:
    :param fragment_file:
    :param apo_file:
    :param dir:
    :return:
    """
    from pymol import cmd
    new_filename = os.path.join(dir, f"{new_fname}.pdb")
    if not os.path.exists(new_filename):
        cmd.reinitialize()
        cmd.load(apo_file, "prot")
        cmd.load(fragment_file, "ligand")
        cmd.create("complex", "ligand, prot")
        cmd.save(new_filename, "complex")
    return new_filename


def calc_interactions(complex_file, interactions_file=None, write_interactions_file=True, return_individual_interactions=False, no_hydro=False):
    """
    Interaction calculation using PLIP (not used)

    :param complex_file:
    :param interactions_file:
    :param write_interactions_file:
    :param return_individual_interactions:
    :param no_hydro:
    :return:
    """
    from plip.basic import config as plip_config
    from plip.structure.preparation import PDBComplex
    if no_hydro:
        plip_config.NOHYDRO = True
        plip_config.VERBOSE = False
        plip_config.QUIET = True
    # if interactions_file and os.path.exists(interactions_file):
    #     interaction_dict = load_json(interactions_file)
    # else:
    mol = PDBComplex()
    mol.load_pdb(complex_file)
    mol.analyze()
    bsid = [":".join([x.hetid, x.chain, str(x.position)]) for x in mol.ligands][0]  # get binding site id
    interactions = mol.interaction_sets[bsid]
    hydrophobic_contacts = [
        f"{int_type.restype}-{int_type.resnr}-0"
        for int_type in interactions.hydrophobic_contacts
    ]
    pi_stacking = [
        f"{int_type.restype}-{int_type.resnr}-1"
        for int_type in interactions.pistacking
    ]
    pi_cation = [
        f"{int_type.restype}-{int_type.resnr}-2"
        for int_type in interactions.pication_laro
    ] + [
        f"{int_type.restype}-{int_type.resnr}-3"
        for int_type in interactions.pication_paro
    ]
    hbond_pdon = [
        f"{int_type.restype}-{int_type.resnr}-4"
        for int_type in interactions.hbonds_pdon
    ]
    hbond_ldon = [
        f"{int_type.restype}-{int_type.resnr}-5"
        for int_type in interactions.hbonds_ldon
    ]
    saltbridge_lneg = [
        f"{int_type.restype}-{int_type.resnr}-6"
        for int_type in interactions.saltbridge_lneg
    ]
    saltbridge_pneg = [
        f"{int_type.restype}-{int_type.resnr}-7"
        for int_type in interactions.saltbridge_pneg
    ]
    saltbridge_metal = [
        f"{int_type.restype}-{int_type.resnr}-8"
        for int_type in interactions.metal_complexes
    ]
    halogen = [
        f"{int_type.restype}-{int_type.resnr}-9"
        for int_type in interactions.halogen_bonds
    ]

    interaction_dict = {
        "hydrophobic_contacts": hydrophobic_contacts,
        "pi_stacking": pi_stacking,
        "pi_cation": pi_cation,
        "hbond_pdon": hbond_pdon,
        "hbond_ldon": hbond_ldon,
        "saltbridge_lneg": saltbridge_lneg,
        "saltbridge_pneg": saltbridge_pneg,
        "saltbridge_metal": saltbridge_metal,
        "halogen": halogen,
    }

    if write_interactions_file:
        dump_json(interaction_dict, interactions_file)

    if return_individual_interactions:
        return hydrophobic_contacts + pi_stacking + pi_cation + hbond_pdon + hbond_ldon + saltbridge_lneg + saltbridge_pneg + saltbridge_metal + halogen

    n_interactions = sum([len(interaction_dict[key]) for key in interaction_dict])
    return n_interactions


def order_best_sucos(mols, fA_mols, fB_mols, sucos_threshold=0.55):
    """
    Order mols, fragment A mols, fragment B mols according to the sucos value (saved as prop) of the mol

    :param mols:
    :param fA_mols:
    :param fB_mols:
    :param sucos_threshold:
    :return:
    """
    # filter all using threshold
    fA_mols = [fA_mol for fA_mol, mol in zip(fA_mols, mols) if float(mol.GetProp('sucos')) >= sucos_threshold]
    fB_mols = [fB_mol for fB_mol, mol in zip(fB_mols, mols) if float(mol.GetProp('sucos')) >= sucos_threshold]
    mols = [mol for mol in mols if float(mol.GetProp('sucos')) >= sucos_threshold]
    # get sucos values
    sucos = [float(mol.GetProp('sucos')) for mol in mols]
    # order all according to the sucos value
    mols = [x for (y, x) in sorted(zip(sucos, mols), reverse=True, key=lambda pair: pair[0])]
    fA_mols = [x for (y, x) in sorted(zip(sucos, fA_mols), reverse=True, key=lambda pair: pair[0])]
    fB_mols = [x for (y, x) in sorted(zip(sucos, fB_mols), reverse=True, key=lambda pair: pair[0])]
    return mols, fA_mols, fB_mols


def get_attach_atom(mol, mcs):
    """
    Get the attachment atom using a molecule and an MCS (by checking if MCS atoms have neighbours that are not part
    of the MCS)

    :param mol:
    :param mcs:
    :return: index of attachment atom
    """
    mcs_atoms = mol.GetSubstructMatch(mcs)
    non_mcs_atoms = [i for i in range(mol.GetNumAtoms()) if i not in mcs_atoms]
    attach_atoms = set()
    for atom in mcs_atoms:
        neighs = [i.GetIdx() for i in mol.GetAtomWithIdx(atom).GetNeighbors()]
        for neigh in neighs:
            if neigh in non_mcs_atoms:
                attach_atoms.add(atom)
    return list(attach_atoms)[0]


def run_attach_atom_check(original_mol, intermed_mol, modified_mol):
    """
    Check the attachment atom in the new molecule is the attachment atom in the intermediate mol

    :param original_mol:
    :param intermed_mol:
    :param modified_mol:
    :return:
    """
    mcs = intermed_mol
    # get attachment atoms in the original molecule and the modified molecule
    original_attach = get_attach_atom(original_mol, mcs)
    new_attach = get_attach_atom(modified_mol, mcs)
    # mcs ordered in same way
    og_matches = original_mol.GetSubstructMatch(mcs)
    new_matches = modified_mol.GetSubstructMatch(mcs)
    # check the attachment atom according to intermed mol same as modified mol
    # (so sub removed and added at same place)
    idx = og_matches.index(original_attach)
    if new_matches[idx] != new_attach:
        return False
    else:
        return True


def run_attach_atom_check_for_loose_expansion(intermed2, intermed1, expansion):
    """
    Check that substructure added to atom in the new linker (difference in atoms between intermed2 and intermed1)
    original_merge -(contraction)-> intermed2 -(expansion)-> intermed1 -(expansion)-> new_merge

    :param intermed2:
    :param intermed1:
    :param expansion:
    :return:
    """
    # get atoms idxs of the new linker atoms
    common_atoms = intermed1.GetSubstructMatch(intermed2)
    new_linker_atoms = [i for i in range(intermed1.GetNumAtoms()) if i not in common_atoms]
    # check that the final expansion is added to one of the new linker atoms
    attach_atom = get_attach_atom(expansion, intermed1)
    if attach_atom in new_linker_atoms:
        return True
    else:
        return False


def calc_prolif_interaction(lig_file, prot_file, lig_protonated=False, interactions=config.PROLIF_INTERACTIONS,
                            interactions_file=None, write_interactions_file=True, ligAsMol=False):
    """
    Calculate interaction using PROLIF, make sure to have hydrogens present

    :param lig_file:
    :param prot_file:
    :param lig_protonated:
    :param interactions:
    :param interactions_file:
    :param write_interactions_file:
    :param ligAsMol:
    :return: list of interactions made
    """
    import prolif as plf
    if interactions_file:
        if os.path.exists(interactions_file):
            try:
                print('ProLIF interactions already calculated in file', interactions_file)
                return load_json(interactions_file)
            except:
                print('faulty file, redoing:', interactions_file)

    if lig_protonated:
        if ligAsMol:
            lig = lig_file
        else:
            lig = Chem.MolFromMolFile(lig_file, removeHs=False)
    else:
        if ligAsMol:
            lig = lig_file
        else:
            lig = Chem.MolFromMolFile(lig_file)
        lig = Chem.AddHs(lig, addCoords=True)
    lig = plf.Molecule.from_rdkit(lig, 'lig')

    try:
        prot = Chem.MolFromPDBFile(prot_file, removeHs=False)
        prot = plf.Molecule.from_rdkit(prot, 'prot')
    except Exception as e:
        print('sanitization error, trying sanitize=False', prot_file)
        prot = Chem.MolFromPDBFile(prot_file, sanitize=False, removeHs=False)
        prot = plf.Molecule.from_rdkit(prot, 'prot')

    fp = plf.Fingerprint(interactions=interactions)
    ifp = fp.generate(lig, prot, metadata=True)
    residue_pairs = list(ifp.keys())

    interaction_list = []
    for residue_pair in residue_pairs:
        residue = f"{residue_pair[1].name}-{residue_pair[1].number}"
        int_types = list(ifp[residue_pair].keys())
        for int_type in int_types:
            for interaction in ifp[residue_pair][int_type]:
                interaction_list.append(f"{int_type}_{residue}")
                # ligand_indices = interaction['parent_indices']['ligand']

    if write_interactions_file:
        dump_json(interaction_list, interactions_file)

    return interaction_list
