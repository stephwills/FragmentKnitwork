from FragmentKnitwork.utils.quilterUtils import check_for_overlap


def check_overlap(mol, substructsA, substructsB):
    """

    :param mol:
    :param substructsA:
    :param substructsB:
    :return:
    """
    overlap_check = False
    for substructA in substructsA:
        for substructB in substructsB:
            overlap_check = check_for_overlap(mol, substructA, substructB)
            if overlap_check:
                break
    return overlap_check
