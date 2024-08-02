import os

from FragmentKnitwork.Quilter.pharmacophore import \
    embed_substructure_with_pharmacophores
from FragmentKnitwork.utils.quilterUtils import (calc_prolif_interaction,
                                                 get_ref_substructures)
from FragmentKnitwork.utils.utils import (disable_rdlogger, dump_json,
                                          get_intersect, get_protein,
                                          order_by_lst)
from joblib import Parallel, delayed
from rdkit import Chem
from tqdm import tqdm

disable_rdlogger()


def _get_ifps(ref_synthon, used_synthon, fB, target, substructure_dir, tmp_dir):
    """
    Get how many interactions are potentially replicated by the replacement substructure with respect to the original
    substructure. Have to iterate through all possible substructure confs (may be multiple matches with the fragment).

    :param ref_synthon:
    :param used_synthon:
    :param fB:
    :param target:
    :param substructure_dir:
    :param tmp_dir:
    :return:
    """
    ref_synthon_files = get_ref_substructures(ref_synthon, fB, isSynthon=True, asMols=False, substructure_dir=substructure_dir)
    prot_file = get_protein(target, fB, protonated=True, desolv=True)

    intersect_ints = []
    for i, ref_synthon_file in enumerate(ref_synthon_files):
        ref_synthon_mol = Chem.MolFromMolFile(ref_synthon_file)
        used_synthon_mol = Chem.MolFromSmiles(used_synthon)
        embed, _, _, _ = embed_substructure_with_pharmacophores(ref_synthon_mol, used_synthon_mol, scoring_mode='max')
        tmp_embed_file = os.path.join(tmp_dir, f"{fB}_{used_synthon}_{ref_synthon}_{i}.mol")
        if embed:
            Chem.MolToMolFile(embed, tmp_embed_file)
            used_ifp = calc_prolif_interaction(tmp_embed_file, prot_file, lig_protonated=False, write_interactions_file=False)
            ref_ifp = calc_prolif_interaction(ref_synthon_file, prot_file, lig_protonated=False, write_interactions_file=False)
            intersect_ints.append(len(get_intersect(used_ifp, ref_ifp)))
        else:
            intersect_ints.append(0)

    if sum(intersect_ints) == 0:
        return None
    else:
        return max(intersect_ints)


def prioritize_data(data, fB, target, substructure_dir, working_dir, n_cpus, new_file, max_run=None):
    """
    Prioritize query results using ProLIF -> we prioritize by merges whereby the replacement substructure is predicted
    to make the same interactions as the substructure it was derived from. This is done by quick ph4 embedding of the
    new substructure using the old, calculate number of overlapping interactions, ruling out those that replicate none
    and ordering by those that make a certain amount.

    :param data:
    :param fB:
    :param target:
    :param substructure_dir:
    :param working_dir:
    :param n_cpus:
    :param new_file:
    :param max_run:
    :return:
    """
    sub_pairs, ref_synthons, used_synthons, expansions, cmpd_ids, similarities = [], [], [], [], [], []

    for substructure_pair in data:
        ref_synthon = substructure_pair.split('*')[1]
        n = len(data[substructure_pair]['expansions'])
        used_synthons.extend(data[substructure_pair]['used_synthons'])
        expansions.extend(data[substructure_pair]['expansions'])
        cmpd_ids.extend(data[substructure_pair]['cmpd_ids'])
        similarities.extend(data[substructure_pair]['similarities'])
        ref_synthons.extend([ref_synthon] * n)
        sub_pairs.extend([substructure_pair] * n)

    # get all the unique inspirational synthons and their proposed replacements
    uniq_synthon_pairs = [f"{ref_synthon}*{used_synthon}" for ref_synthon, used_synthon in
                          zip(ref_synthons, used_synthons)]
    uniq_synthon_pairs = list(set(uniq_synthon_pairs))
    uniq_synthon_pairs.sort()
    uniq_refs = [p.split('*')[0] for p in uniq_synthon_pairs]
    uniq_used = [p.split('*')[1] for p in uniq_synthon_pairs]

    # calculate their IFPs to to allow prioritization by which retains interactions the best
    print('Calculating IFPs')
    res = Parallel(n_jobs=n_cpus, backend="multiprocessing")(
        delayed(_get_ifps)(
            ref_synthon, used_synthon, fB, target, substructure_dir, working_dir
        ) for ref_synthon, used_synthon in tqdm(zip(uniq_refs, uniq_used), total=len(uniq_refs))
    )
    print('IFPs calculated, running processing')

    uniq_int_counts = {pair: int_count for pair, int_count in zip(uniq_synthon_pairs, res)}
    int_counts = []
    for ref_synthon, used_synthon in zip(ref_synthons, used_synthons):
        P = f"{ref_synthon}*{used_synthon}"
        int_counts.append(uniq_int_counts[P])

    def filt_by_bool(to_filter, bool_list):
        return [i for i, b in zip(to_filter, bool_list) if b]

    sub_pairs = filt_by_bool(sub_pairs, int_counts)
    used_synthons = filt_by_bool(used_synthons, int_counts)
    expansions = filt_by_bool(expansions, int_counts)
    cmpd_ids = filt_by_bool(cmpd_ids, int_counts)
    similarities = filt_by_bool(similarities, int_counts)
    int_counts = [i for i in int_counts if i]

    if max_run:
        N = max_run
    else:
        N = len(expansions)

    sub_pairs = order_by_lst(sub_pairs, int_counts)[:N]
    used_synthons = order_by_lst(used_synthons, int_counts)[:N]
    expansions = order_by_lst(expansions, int_counts)[:N]
    cmpd_ids = order_by_lst(cmpd_ids, int_counts)[:N]
    similarities = order_by_lst(similarities, int_counts)[:N]

    uniq_sub_pairs = list(data.keys())
    new_data = {pair: {'expansions': [],
                       'used_synthons': [],
                       'cmpd_ids': [],
                       'similarities': []} for pair in uniq_sub_pairs}
    for sub_pair, used_synthon, exp, cmpd_id, sim in zip(sub_pairs, used_synthons, expansions, cmpd_ids, similarities):
        new_data[sub_pair]['expansions'].append(exp)
        new_data[sub_pair]['used_synthons'].append(used_synthon)
        new_data[sub_pair]['cmpd_ids'].append(cmpd_id)
        new_data[sub_pair]['similarities'].append(sim)

    dump_json(new_data, new_file)
    return new_data
