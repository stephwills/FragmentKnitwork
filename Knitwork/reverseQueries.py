import os
from argparse import ArgumentParser
import shutil
import time

from FragmentKnitwork.Knitwork.queries import replace_substructure
from FragmentKnitwork.utils import knitworkConfig as config
from FragmentKnitwork.utils.knitworkUtils import get_property_dict
from FragmentKnitwork.utils.utils import load_json, dump_json
from joblib import Parallel, delayed
from rdkit import Chem
from tqdm import tqdm


def find_analogues(name, subnode_equivalent_data, property_dict, working_dir, threshold=config.SIMILARITY_THRESHOLD, limit=config.REVERSE_QUERY_LIMIT, strict=config.REVERSE_QUERY_STRICT):
    """
    Run the reverse queries

    :param name: name of the pair
    :param subnode_equivalent_data: data to get the equivalent synthon pattern for a subnode (with attachment point specified)
    :param property_dict: property dict for mols (because it isn't preserved when passing around)
    :param working_dir: where to save files
    :param threshold: similarity threshold
    :param limit: limit number of results
    :param strict: whether to allow change in linker
    :return: expansion data
    """
    expansions_dict = {}

    output_json = os.path.join(working_dir, f"{name}_expansions.json")
    if strict:
        timings_file = os.path.join(working_dir, f"{name}_reverse_merge_strict_timings.json")
    else:
        timings_file = os.path.join(working_dir, f"{name}_reverse_merge_timings.json")
    if os.path.exists(output_json):
        print('Already run for mol', name)
        expansions_dict = load_json(output_json)
        return expansions_dict

    start = time.time()

    # load property names from dict (because they are lost when passing mol between functions)
    fragmentA = property_dict[name]['fragmentA']
    subnode = property_dict[name]['ref_subnode']
    synthon = property_dict[name]['ref_synthon']
    used_synthon = property_dict[name]['used_synthon']
    original_similarity = property_dict[name]['similarity']
    original_smiles = property_dict[name]['smiles']

    # get equivalent synthon # TODO: probably better way to do this in future
    if subnode not in subnode_equivalent_data[fragmentA]:
        print('Could not find subnode', subnode, 'in equivalent subnode file')
        return None

    # get the SMILES for the initial substructure we want to replace IN form where the equivalent edge label is used
    # the edge labels have attachment points unlike SMILES on nodes
    ref_synthons = subnode_equivalent_data[fragmentA][subnode].split(',')

    # there may be multiple possible synthons
    num_expansions = 0
    for i, ref_synthon in enumerate(ref_synthons):
        key = f"{ref_synthon}*{synthon}*{original_smiles}"
        # run the reverse query
        rev_expansions, rev_synthons, rev_similarities, rev_cmpd_ids, rev_intermed_nodes1, rev_intermed_nodes2 = replace_substructure(
            original_smiles, ref_synthon, threshold, limit, strict, used_synthon)

        if len(rev_expansions) > 0:
            num_expansions += len(rev_expansions)
            inner_dict = {'expansions': rev_expansions,
                          'cmpd_ids': rev_cmpd_ids,
                          'names': [f"{name}-{idx}" for idx in range(len(rev_expansions))],
                          'used_synthons': [used_synthon] * len(rev_expansions),
                          'rev_synthons': rev_synthons,
                          'similarities': rev_similarities,
                          'intermed_nodes1': rev_intermed_nodes1,
                          'original_names': [name] * len(rev_expansions),
                          'original_similarities': [float(original_similarity)] * len(rev_expansions)}
            if not strict:
                inner_dict['intermed_nodes2'] = rev_intermed_nodes2
            expansions_dict[key] = inner_dict

    dump_json(expansions_dict, output_json)
    if num_expansions != 0:
        print(f'Finished querying for smiles: {original_smiles}, name {name}, num expansions:', num_expansions)

    end = time.time()
    time_taken = round(end-start, 2)
    dump_json(time_taken, timings_file)

    return expansions_dict


def runReverseMerging(sdf_file, subnode_equivalent_file, max_run, n_parallel, output_dir, working_dir, threshold, limit_results, strict_linker_search):
    """
    Performs another round of bioisosteric/impure merging of already placed impure merges. Provide with the SDF file containing scored bioisosteric merges.

    :param sdf_file: sdf file containing scored impure merges
    :param subnode_equivalent_file: generated during runEnumeration.py
    :param max_run: max number of molecules to run reverse merging on
    :param n_parallel: number of parallel queries
    :param output_dir: output dir
    :param working_dir: where to save intermediate files
    :param threshold: similarity threshold for pharm similarity calculation
    :param limit_results: limit number of possible reverse merges for each first round bioisosteric merge
    :param strict_linker_search: whether to allow flexibility in the linker or not
    :return:
    """
    subnode_equivalent_data = load_json(subnode_equivalent_file)

    # load all mols from the sdf file
    mols = list(Chem.SDMolSupplier(sdf_file))

    # optional limit how many molecules run
    if not max_run:
        max_run = len(mols)
    else:
        max_run = max_run
    mols = mols[:max_run]

    # each fragment pair has a file
    pairs = [f"{mol.GetProp('fragmentA')}-{mol.GetProp('fragmentB')}" for mol in mols]
    names = [mol.GetProp('_Name') for mol in mols]
    uniq_pairs = list(set(pairs))

    # create a dictionary for each pair {pair: [mols]}
    mol_dict = {pair: [] for pair in uniq_pairs}
    for name, pair in zip(names, pairs):
        mol_dict[pair] += [name]

    # for some reason the properties aren't preserved when running parallelisation
    property_dict = get_property_dict(mols)  # {name: {prop_name: prop, prop_name: prop}}

    start = time.time()
    # run the queries
    results = Parallel(n_jobs=n_parallel, backend="multiprocessing")(
        delayed(find_analogues)(name, subnode_equivalent_data, property_dict, working_dir, threshold, limit_results, strict_linker_search)
        for name in tqdm(names, position=0, leave=True, total=len(names))
    )
    print('Queries run')
    end = time.time()

    total_time = round(end-start, 2)

    # process results into individual files for fragment pairs
    for uniq_pair in uniq_pairs:
        print(uniq_pair)
        pair_file = os.path.join(output_dir, f"{uniq_pair}_reverse_merge.json")
        if not os.path.exists(pair_file):
            pair_expansions_dict = {}
            for pair, res_dict in zip(pairs, results):
                if uniq_pair == pair:
                    if res_dict:
                        if len(res_dict) > 0:
                            pair_expansions_dict.update(res_dict)

            dump_json(pair_expansions_dict, pair_file)
        else:
            print(pair_file, 'already exists -- please check!')

    time_data = {'total_time': total_time,
                 'n_cpus': n_parallel}
    if strict_linker_search:
        dump_json(time_data, os.path.join(output_dir, 'reverse_merge_strict_total_time.json'))
    else:
        dump_json(time_data, os.path.join(output_dir, 'reverse_merge_loose_total_time.json'))


def main():
    """
    This is code to further diversity the set of bioisosteric merges by taking the merges and replacing the other
    substructure used to form the merge (represented by the initial seed node in the substructure pair for the
    impure_merging queries). It can either be run strict (no change in linker) or allow some change in the linker. We
    use a processed filtered SDF as an input from the placed impure merges (to focus on the most promising).

    :return:
    """
    parser = ArgumentParser()
    parser.add_argument('--sdf_file')
    parser.add_argument('--substructure_dir', help='dir where all the substructure data has been saved generated in runEnumeration.py')
    parser.add_argument('--working_dir', help='to save intermediate files')
    parser.add_argument('--output_dir', help='to save output files')
    parser.add_argument('--n_parallel', type=int)
    parser.add_argument('--limit_results', type=int, default=config.REVERSE_QUERY_LIMIT)
    parser.add_argument('--threshold', type=float, default=config.SIMILARITY_THRESHOLD)
    parser.add_argument('--linker_option', choices=['strict', 'loose', 'both'], help='whether to allow no change in linker, change in linker, or run both options')
    parser.add_argument('--max_run', type=int, default=None, help='if we dont want to run for every molecule')
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    if not os.path.exists(args.working_dir):
        os.mkdir(args.working_dir)

    subnode_equivalent_file = os.path.join(args.substructure_dir, 'equivalent_subnodes.json')
    if args.linker_option == 'strict':
        runReverseMerging(args.sdf_file, subnode_equivalent_file, args.max_run, args.n_parallel, args.output_dir, args.working_dir,
                          args.threshold, args.limit_results, True)
    if args.linker_option == 'loose':
        runReverseMerging(args.sdf_file, subnode_equivalent_file, args.max_run, args.n_parallel, args.output_dir, args.working_dir,
                          args.threshold, args.limit_results, False)
    if args.linker_option == 'both':
        strict_work_dir, loose_work_dir = os.path.join(args.working_dir, 'strict_search'), os.path.join(args.working_dir, 'loose_search')
        strict_out_dir, loose_out_dir = os.path.join(args.ouptut_dir, 'strict_search'), os.path.join(args.output_dir, 'loose_search')
        os.mkdir(strict_work_dir)
        os.mkdir(strict_out_dir)
        os.mkdir(loose_work_dir)
        os.mkdir(loose_out_dir)
        print('Running strict linker search')
        runReverseMerging(args.sdf_file, subnode_equivalent_file, args.max_run, args.n_parallel, strict_out_dir, strict_work_dir,
                          args.threshold, args.limit_results, True)
        print('Running loose linker search')
        runReverseMerging(args.sdf_file, subnode_equivalent_file, args.max_run, args.n_parallel, loose_out_dir, loose_work_dir,
                          args.threshold, args.limit_results, False)

    # shutil.rmtree(args.working_dir)


if __name__ == "__main__":
    main()
