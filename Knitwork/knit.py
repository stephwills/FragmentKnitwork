"""Code for running by fragment pair rather than by substructure pair (not actively used)"""
import os

from FragmentKnitwork.Knitwork.queries import impure_expansion, pure_expansion
from FragmentKnitwork.utils.knitworkUtils import calc_pharm_fp, load_sigFactory
from FragmentKnitwork.utils.utils import dump_json, load_json
from rdkit import Chem

######## THIS CODE NOT ACTIVELY USED - FOR RUNNING BY FRAGMENT PAIR RATHER THAN SUBSTRUCTURE ########
def pure_merge(fragment_pair, substructure_pairs, output_dir, write_to_file=True, return_results=False, limit=None):
    """

    :param fragment_pair:
    :param substructure_pairs:
    :param output_dir:
    :param write_to_file:
    :param return_results:
    :return:
    """
    print('STARTING fragment pair', fragment_pair)
    output_fname = os.path.join(output_dir, f"{fragment_pair}_pure_merge.json")
    if os.path.exists(output_fname):
        print('File already exists for these queries.')
        if return_results:
            return load_json(output_fname)
        else:
            return None

    output_data = {}

    for i, substructure_pair in enumerate(substructure_pairs):
        print(f"{fragment_pair}:", 'Running substructure pair', substructure_pair, f"{i + 1}/{len(substructure_pairs)}")
        subnode, synthon = substructure_pair[0], substructure_pair[1]

        expansions, cmpd_ids = pure_expansion(subnode, synthon, limit=limit)
        print(f"{fragment_pair}:", 'Found', len(expansions), 'expansions')
        substructure_pair_data = {'expansions': expansions,
                                  'cmpd_ids': cmpd_ids}
        key = f"{subnode}*{synthon}"
        output_data[key] = substructure_pair_data

    print('FINISHED fragment pair', fragment_pair)
    if write_to_file:
        dump_json(output_data, output_fname)

    if return_results:
        return output_data


def impure_merge(fragment_pair, substructure_pairs, desc, output_dir, write_to_file=True, return_results=False, limit=None):
    """

    :param fragment_pair:
    :param substructure_pairs:
    :param desc:
    :param output_dir:
    :param write_to_file:
    :param return_results:
    :return:
    """
    print('STARTING fragment pair', fragment_pair)
    output_fname = os.path.join(output_dir, f"{fragment_pair}_{desc}_impure_merge.json")
    if os.path.exists(output_fname):
        print('File already exists for these queries.')
        if return_results:
            return load_json(output_fname)
        else:
            return None

    if desc == 'prop_pharmfp':
        sigFactory = load_sigFactory()

    output_data = {}

    for i, substructure_pair in enumerate(substructure_pairs):
        print(f"{fragment_pair}:", 'Running substructure pair', substructure_pair, f"{i+1}/{len(substructure_pairs)}")
        subnode, synthon = substructure_pair[0], substructure_pair[1]
        if desc == 'prop_pharmfp':
            vector = calc_pharm_fp(Chem.MolFromSmiles(synthon), sigFactory, asStr=False)
        expansions, repl_synthons, similarities, cmpd_ids = impure_expansion(subnode, vector, synthon, limit=limit)
        print(f"{fragment_pair}:", 'Found', len(expansions), 'expansions')
        substructure_pair_data = {'expansions': expansions,
                                  'used_synthons': repl_synthons,
                                  'similarities': similarities,
                                  'cmpd_ids': cmpd_ids}

        key = f"{subnode}*{synthon}"
        output_data[key] = substructure_pair_data

    print('FINISHED fragment pair', fragment_pair)
    if write_to_file:
        dump_json(output_data, output_fname)

    if return_results:
        return output_data
