from FragmentKnitwork.utils.quilterUtils import calc_prolif_interaction
from FragmentKnitwork.utils.utils import get_intersect, split_complex_file


def check_interactions(holo_file, mol_file, fA, fB, fragment_ints_data, interactions_file):
    """

    :param holo_file:
    :param mol_file:
    :param fA:
    :param fB:
    :param fragment_ints_data:
    :param interactions_file:
    :return:
    """
    # get the pre-calculated fragment interactions
    fA_ints, fB_ints = fragment_ints_data[fA], fragment_ints_data[fB]
    uniq_fA_ints, uniq_fB_ints = list(set(fA_ints)), list(set(fB_ints))

    # get the unique fragment interactions from both fragments
    uniq_frag_ints = list(set(uniq_fA_ints + uniq_fB_ints))

    # get apo file (for prolif calculation)
    prot_file = holo_file.replace('.holo_minimised.pdb', '.apo_minimised.pdb')
    prot_file = split_complex_file(holo_file, prot_file)

    # calculate the molecule interactions
    mol_ints = calc_prolif_interaction(mol_file, prot_file, lig_protonated=True, interactions_file=interactions_file,
                                       write_interactions_file=True, ligAsMol=False)
    uniq_mol_ints = list(set(mol_ints))

    # save molecule interactions as string
    str_mol_ints = ','.join(mol_ints)

    # total number of interactions made
    total_ints = len(mol_ints)
    # total number of interactions made (not include hydrophobic)
    total_ints_nonhyd = len([i for i in uniq_mol_ints if 'Hydrophobic' not in i])

    # calculate various metrics
    num_uniq_frag_ints = len(uniq_frag_ints)
    # calculate the number of fragment interactions that are replicated (at residue level)
    num_repl_uniq_ints = len(get_intersect(uniq_frag_ints, uniq_mol_ints))

    # calculate number of interactions lost and that are new
    num_uniq_ints_lost = len([i for i in uniq_frag_ints if i not in uniq_mol_ints])
    num_new_uniq_ints = len([i for i in uniq_mol_ints if i not in uniq_frag_ints])

    # check that an interaction seen in fA and in fB are replicated for the main int check
    fA_ints_rep = len(get_intersect(uniq_fA_ints, uniq_mol_ints)) > 0
    fB_ints_rep = len(get_intersect(uniq_fB_ints, uniq_mol_ints)) > 0
    int_check = sum([fA_ints_rep, fB_ints_rep]) == 2

    int_data = {'interactions': str_mol_ints,
                'total_num_ints': total_ints,
                'total_num_ints_nonhyd': total_ints_nonhyd,
                'num_uniq_frag_ints': num_uniq_frag_ints,
                'num_repl_uniq_ints': num_repl_uniq_ints,
                'num_uniq_ints_lost': num_uniq_ints_lost,
                'num_new_uniq_ints': num_new_uniq_ints,
                'int_check': int_check}

    return int_check, int_data
