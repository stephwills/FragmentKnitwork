"""Perform R group expansion of a given merge"""


from FragmentKnitwork.Knitwork.queries import single_expansion
from FragmentKnitwork.utils.utils import dump_json


def get_relevant_r_groups(_ref_subnode, ref_synthon, fragmentA, fragmentB, r_group_data, equiv_synthon_data):
    """
    Gets a bit complicated - get the possible R groups to used in expansion queries

    :param _ref_subnode:
    :param ref_synthon:
    :param fragmentA:
    :param fragmentB:
    :param r_group_data:
    :param equiv_synthon_data:
    :return:
    """
    if '[Xe]' in _ref_subnode:
        ref_subnodes = [_ref_subnode]
    else:
        if _ref_subnode in equiv_synthon_data[fragmentA]:
            ref_subnodes = equiv_synthon_data[fragmentA][_ref_subnode].split(',')
        else:
            print(_ref_subnode, 'not in equiv synthon dict')
            return None

    fA_expansions = []
    fB_expansions = []
    for ref_subnode in ref_subnodes:
        if ref_subnode in r_group_data[fragmentA]:
            _fA_expansions = r_group_data[fragmentA][ref_subnode]
            fA_expansions.extend(_fA_expansions)
    if ref_synthon in r_group_data[fragmentB]:
        fB_expansions = r_group_data[fragmentB][ref_synthon]

    r_groups = fA_expansions + fB_expansions
    r_groups = list(set(r_groups))
    return r_groups


def run_r_group_expansion(sub_pair_data, sub_pair, _ref_subnode, ref_synthon, fragmentA, fragmentB, r_group_data, equiv_synthon_data):
    """
    Run R-group expansions during the querying stage for pure or impure (first round) merges - will result in quite a lot of data!
    Adds to the existing data dictionary and writes to new separate output directory.

    :param merges:
    :param _ref_subnode:
    :param ref_synthon:
    :param fragmentA:
    :param fragmentB:
    :param r_group_data:
    :param equiv_synthon_data:
    :return:
    """
    total_expansions = 0
    new_names = ['expansions', 'cmpd_ids', 'r_group', 'smiles_before_rgrp_expansion']
    property_names = [name for name in list(sub_pair_data.keys()) if name not in new_names]
    merges = sub_pair_data['expansions']
    sub_pair_data['r_group'] = [None] * len(merges)
    sub_pair_data['smiles_before_rgrp_expansion'] = [None] * len(merges)

    r_groups = get_relevant_r_groups(_ref_subnode, ref_synthon, fragmentA, fragmentB, r_group_data, equiv_synthon_data)

    for i, merge in enumerate(merges):
        for r_group in r_groups:
            expansions, cmpd_ids = single_expansion(merge, r_group)
            total_expansions += len(expansions)
            if len(expansions) > 0:
                for expansion, cmpd_id in zip(expansions, cmpd_ids):
                    sub_pair_data['expansions'].append(expansion)
                    sub_pair_data['cmpd_ids'].append(cmpd_id)
                    sub_pair_data['r_group'].append(r_group)
                    sub_pair_data['smiles_before_rgrp_expansion'].append(merge)
                    for prop in property_names:
                        sub_pair_data[prop].append(sub_pair_data[prop][i])

    print(f"{sub_pair}: {total_expansions} R group expansions found")
    return sub_pair_data


def write_r_group_expansions(fragment_data, fragment_pair, r_group_data, equiv_synthon_data, r_group_output_file):
    """

    :param fragment_data:
    :param fragment_pair:
    :param r_group_data:
    :param equiv_synthon_data:
    :return:
    """
    r_group_fragment_data = {}
    for sub_pair in fragment_data:
        sub_pair_data = fragment_data[sub_pair]
        new_sub_pair_data = run_r_group_expansion(sub_pair_data, sub_pair, sub_pair.split('*')[0],
                                                  sub_pair.split('*')[1],
                                                  fragment_pair.split('-')[0], fragment_pair.split('-')[1],
                                                  r_group_data,
                                                  equiv_synthon_data)
        r_group_fragment_data[sub_pair] = new_sub_pair_data
    dump_json(r_group_fragment_data, r_group_output_file)
