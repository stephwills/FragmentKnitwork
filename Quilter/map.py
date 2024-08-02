"""Some very complicated code for generating possible mappings using pharmacophore embedding"""
import FragmentKnitwork.utils.quilterConfig as config
from FragmentKnitwork.Quilter.pharmacophore import (
    check_carbons, embed_substructure_with_pharmacophores, loop_list)
from FragmentKnitwork.utils.quilterUtils import coords_from_ids
from FragmentKnitwork.utils.utils import get_distance, get_intersect
from rdkit import Chem
from rdkit.Chem import rdFMCS


def get_matching_coords(repl_substruct, ref_substruct, repl_atom_ids, ref_atom_ids,
                        atom_dist_threshold=config.ATOM_MATCH_THRESHOLD):
    """
    Get corresponding atoms to use to create the atom mapping between the repl_substruct and ref_substruct.

    :param repl_substruct:
    :param ref_substruct:
    :param repl_atom_ids:
    :param ref_atom_ids:
    :param atom_dist_threshold: the distance threshold to use
    :return:
    """
    atom_matches = {}
    for repl_atoms, ref_atoms in zip(repl_atom_ids, ref_atom_ids):
        # if only one atom, then just map to each other
        if len(repl_atoms) == 1 and len(ref_atoms) == 1:
            atom_matches[ref_atoms[0]] = repl_atoms[0]

        else:
            # get coordinates for the repl and ref atoms
            repl_coords, ref_coords = coords_from_ids(repl_atoms, repl_substruct), coords_from_ids(ref_atoms, ref_substruct)
            repl_matches, ref_matches = [], []
            dists, atoms = [], []

            # calculate distances between all pairs of atoms
            for repl_atom, repl_coord in zip(repl_atoms, repl_coords):
                for ref_atom, ref_coord in zip(ref_atoms, ref_coords):
                    dists.append(get_distance(repl_coord, ref_coord))
                    atoms.append((repl_atom, ref_atom))

            # get the closest pairs of atoms that are also below the distance threshold
            ordered_atoms = [x for (y, x) in sorted(zip(dists, atoms), key=lambda pair: pair[0])]
            dists.sort()

            # only get the atom pairs where the distance is below a certain threshold
            n_ats = min([len(repl_atoms), len(ref_atoms)])
            use_atoms = [atom_pair for atom_pair, dist in zip(ordered_atoms[:n_ats], dists[:n_ats]) if dist <= atom_dist_threshold]

            # don't use overlapping atoms (e.g. atom in one sub close to two atoms in other sub)
            for atom_pair in use_atoms:
                if atom_pair[0] not in repl_matches and atom_pair[1] not in ref_matches:
                    repl_matches.append(atom_pair[0])
                    ref_matches.append(atom_pair[1])

            for ref_match, repl_match in zip(ref_matches, repl_matches):
                atom_matches[ref_match] = repl_match

    return [atom_matches]


def get_atom_matches_with_merge(merge, repl_substruct, sub_atom_maps, isSynthon=False):
    """
    Using the mappings between substructures, transfer into map with the full merge (numbers sadly change)

    :param merge:
    :param repl_substruct:
    :param sub_atom_maps:
    :param isSynthon:
    :return:
    """
    all_merge_atom_maps = []
    for sub_atom_map in sub_atom_maps:  # for every possible substructure map (multiple if dealing with rings)

        # get the equivalent numbering between the substructure and merge (ordered in the same way)
        if isSynthon:  # if there is a synthon, have to get MCS because of the attachment (xenon) atom
            mcs = Chem.MolFromSmarts(rdFMCS.FindMCS([merge, repl_substruct]).smartsString)
            merge_substruct_matches = merge.GetSubstructMatches(mcs)  # may be multiple
            repl_substruct_match = repl_substruct.GetSubstructMatch(mcs)

        else:
            merge_substruct_matches = merge.GetSubstructMatches(repl_substruct)
            repl_substruct_match = repl_substruct.GetSubstructMatch(repl_substruct)

        merge_atom_maps = []

        for merge_substruct_match in merge_substruct_matches:
            merge_atom_map = {}
            # for each atom in the map (with substructure), get the equivalent atom idx in the merge
            for key, val in sub_atom_map.items():
                idx = repl_substruct_match.index(val)
                merge_atom = merge_substruct_match[idx]
                merge_atom_map[key] = merge_atom

            merge_atom_maps.append(merge_atom_map)  # for all possible substructure matches
        all_merge_atom_maps.extend(merge_atom_maps)  # for all possible orderings

    return all_merge_atom_maps


def get_sub_custom_map(ref_sub, used_sub, merge, isSynthon=False):
    """

    :param ref_sub:
    :param used_sub:
    :param merge:
    :param isSynthon:
    :return:
    """
    embedded, ref_ids, used_ids, ring_map = embed_substructure_with_pharmacophores(ref_sub, used_sub)

    if not embedded:
        return None

    if ring_map:
        atom_matches_with_sub = ring_map
    else:
        atom_matches_with_sub = get_matching_coords(embedded, ref_sub, used_ids, ref_ids)

    # atom matches with sub is a LIST of DICTIONARIES giving the atom mapping
    # it is a LIST because with ring map there may be multiple possible mappings
    merge_matches_with_sub = get_atom_matches_with_merge(merge, used_sub, atom_matches_with_sub, isSynthon=isSynthon)
    return merge_matches_with_sub


def get_custom_maps(merge, ref_subnode, ref_synthon, used_subnode, used_synthon, reverse_merge=False):
    """
    Get possible custom maps between the ref substructures used for a merge and the merge.

    :param merge:
    :param ref_subnode:
    :param ref_synthon:
    :param used_subnode:
    :param used_synthon:
    :return:
    """
    merge_matches_with_synthon, merge_matches_with_subnode = None, None
    if not reverse_merge:
        merge_matches_with_subnode = get_sub_custom_map(ref_subnode, used_subnode, merge)
        merge_matches_with_synthon = get_sub_custom_map(ref_synthon, used_synthon, merge, isSynthon=True)
    if reverse_merge:
        merge_matches_with_subnode = get_sub_custom_map(ref_subnode, used_subnode, merge, isSynthon=True)
        merge_matches_with_synthon = get_sub_custom_map(ref_synthon, used_synthon, merge, isSynthon=True)

    if not merge_matches_with_synthon or not merge_matches_with_subnode:
        return None

    custom_maps = []

    for subnode_map in merge_matches_with_subnode:
        for synthon_map in merge_matches_with_synthon:
            # check the atoms aren't overlapping
            if len(get_intersect(list(subnode_map.values()), list(synthon_map.values()))) == 0:
                custom_maps.append((subnode_map, synthon_map))

    return custom_maps


def check_carbon_rings_in_fragment_custom_maps(feats, fragment_ids, fragment, merge_ids, merge):
    """

    :param feats:
    :param fragment_ids:
    :param fragment:
    :param merge_ids:
    :param merge:
    :return:
    """
    rings = ['Aliphatic', 'Aromatic']
    custom_maps = []
    if len(feats) == 1 and feats[0].GetFamily() in rings:
        if check_carbons(fragment_ids[0], fragment) or check_carbons(merge_ids[0], merge):
            looped_fragment_ids = loop_list(fragment_ids[0])
            for loop in looped_fragment_ids:
                d = {frag_id: merge_id for frag_id, merge_id in zip(loop, merge_ids[0])}
                custom_maps.append(d)
        else:
            d = {frag_id: merge_id for frag_id, merge_id in zip(fragment_ids[0], merge_ids[0])}
            custom_maps.append(d)
    else:
        d = {}
        for _frag_ids, _merge_ids in zip(fragment_ids, merge_ids):
            for frag_id, merge_id in zip(_frag_ids, _merge_ids):
                d[frag_id] = merge_id
        custom_maps.append(d)
    return custom_maps


def get_custom_maps_for_fragments(fragmentA, fragmentB, fragment_atom_ids, merge_atom_ids, feats, merge):
    """

    :param fragmentA:
    :param fragmentB:
    :param fragment_atom_ids:
    :param merge_atom_ids:
    :param feats:
    :param merge:
    :return:
    """
    fA_feats = []
    fA_fragment_ids = []
    fA_merge_ids = []

    fB_feats = []
    fB_fragment_ids = []
    fB_merge_ids = []

    for _fragment_atom_ids, _merge_atom_ids, feat in zip(fragment_atom_ids, merge_atom_ids, feats):
        if feat.GetMol().GetProp('_Name') == 'fA':
            fA_feats.append(feat)
            fA_fragment_ids.append(_fragment_atom_ids)
            fA_merge_ids.append(_merge_atom_ids)

        if feat.GetMol().GetProp('_Name') == 'fB':
            fB_feats.append(feat)
            fB_fragment_ids.append(_fragment_atom_ids)
            fB_merge_ids.append(_merge_atom_ids)

    custom_maps_A = check_carbon_rings_in_fragment_custom_maps(fA_feats, fA_fragment_ids, fragmentA, fA_merge_ids, merge)
    custom_maps_B = check_carbon_rings_in_fragment_custom_maps(fB_feats, fB_fragment_ids, fragmentB, fB_merge_ids, merge)

    maps = []
    for mapA in custom_maps_A:
        for mapB in custom_maps_B:
            maps.append({'fA': mapA,
                         'fB': mapB})

    return maps
