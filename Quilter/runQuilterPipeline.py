
from argparse import ArgumentParser
import os

from FragmentKnitwork.Quilter.quilterPipeline import Pipeline
from FragmentKnitwork.Quilter.quilterWictorPipeline import WictorPipeline
from FragmentKnitwork.utils import quilterConfig as config
from FragmentKnitwork.utils.utils import load_json


def main():
    parser = ArgumentParser()
    parser.add_argument('--json_file', required=True, help='This should be one of the json files generated for a fragment pair by runKnitting.py')
    parser.add_argument('--pdb_file', required=True, help='apo-desolv pdb file to use for placing the merges')
    parser.add_argument('--fragmentA', help='The name of the first fragment in the pair (how it is named in the Fragalysis format)')
    parser.add_argument('--fragmentB', help='The name of the second fragment in the pair (how it is named in the Fragalysis format)')
    parser.add_argument('--target', required=False, help='The name of the target (how it is named in the Fragalysis format)')
    parser.add_argument('--substructure_dir', required=False, default=config.SUBSTRUCTURE_DIR, help='The dir generated by runEnumeration.py')
    parser.add_argument('--output_dir', required=False, default=config.OUTPUT_DIR, help='where to save output data')
    parser.add_argument('--working_dir', required=False, default=config.WORKING_DIR, help='where to save intermediate files')
    parser.add_argument('--scoring_dir', required=False, help='where to save the scored ligands')
    parser.add_argument('--n_cpus', type=int, required=False, default=config.N_CPUS, help='the number of CPUs to parallelise over')
    parser.add_argument('--type_merge', choices=['pure_merge', 'impure_merge', 'reverse_merge'],
                        help='The type of merge, a pure (perfect) merge, or a bioisosteric merge where one (impure merge) or two substructures (reverse merge) has been replaced')
    parser.add_argument('--limit_num_run', action='store_true', default=config.LIMIT_NUM_RUN, help='whether to limit the number of compounds overall filtered')
    parser.add_argument('--limit_num_minimize', action='store_true',
                        help='only relevant if using wictor pre-filter, whether to limit the number of compounds minimized (after the first pre-filter)')
    parser.add_argument('--max_num_filter', type=int, required=False, default=config.MAX_NUM_FILTER, help='the max number of compounds to filter overall')
    parser.add_argument('--max_num_minimize', type=int, required=False, help='only relevant if using wictor pre-filter, the max number of compounds to minimize')
    parser.add_argument('--move_files', action='store_true', default=config.MOVE_FILES, help='(recommended) whether to move files back to the output directory')
    parser.add_argument('--use_wictor', action='store_true', help='whether to use the wictor pre-filter (without minimization) then prioritize the top X')
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    if not os.path.exists(args.working_dir):
        os.mkdir(args.working_dir)

    if not os.path.exists(args.scoring_dir):
        os.mkdir(args.scoring_dir)

    # read in input data
    file = args.json_file
    data = load_json(file)
    smiles, cmpd_ids = [], []
    ref_subnodes, ref_synthons = [], []

    # only used for impure merges & reverse merges
    similarities, used_subnodes, used_synthons = [], [], []

    # only used for reverse merges
    original_smiles, original_names, original_similarities = [], [], []
    rev_synthons, intermed_nodes1, intermed_nodes2 = [], [], []

    # read data from input file
    for sub_pair in data:
        print('Substructure pair', sub_pair)
        n_expans = len(data[sub_pair]['expansions'])
        print('N expansions', n_expans)

        # shared between all merge types
        smiles.extend(data[sub_pair]['expansions'])
        cmpd_ids.extend(data[sub_pair]['cmpd_ids'])

        # get ref substructures (not for reverse merge)
        subnode, synthon = sub_pair.split('*')[0], sub_pair.split('*')[1]
        ref_subnodes.extend([subnode] * n_expans)
        ref_synthons.extend([synthon] * n_expans)

        # if impure merges
        if args.type_merge == 'impure_merge' or args.type_merge == 'reverse_merge':
            similarities.extend(data[sub_pair]['similarities'])
            used_subnodes.extend([subnode] * n_expans)
            used_synthons.extend(data[sub_pair]['used_synthons'])

        if args.type_merge == 'reverse_merge':
            original_smi = sub_pair.split('*')[2]
            original_smiles.extend([original_smi] * n_expans)
            original_similarities.extend(data[sub_pair]['original_similarities'])
            original_names.extend(data[sub_pair]['original_names'])
            rev_synthons.extend(data[sub_pair]['rev_synthons'])
            intermed_nodes1.extend(data[sub_pair]['intermed_nodes1'])
            if 'intermed_nodes2' in list(data[sub_pair].keys()):
                intermed_nodes2.extend(data[sub_pair]['intermed_nodes2'])

    if len(intermed_nodes2) == 0:
        intermed_nodes2 = None

    if args.limit_num_run:
        limit = args.max_num_filter
    else:
        limit = None

    if args.limit_num_minimize:
        minimize_limit = args.max_num_minimize
    else:
        minimize_limit = None

    if not args.use_wictor:
        pipeline = Pipeline(smiles,
                            args.fragmentA,
                            args.fragmentB,
                            args.target,
                            args.pdb_file,
                            cmpd_ids,
                            ref_subnodes,
                            ref_synthons,
                            args.output_dir,
                            args.working_dir,
                            args.scoring_dir,
                            args.type_merge,
                            args.substructure_dir,
                            similarities,
                            used_subnodes,
                            used_synthons,
                            original_similarities,
                            original_names,
                            rev_synthons,
                            intermed_nodes1,
                            intermed_nodes2,
                            limit,
                            args.n_cpus,
                            args.move_files)
    else:
        pipeline = WictorPipeline(smiles,
                                  args.fragmentA,
                                  args.fragmentB,
                                  args.target,
                                  args.pdb_file,
                                  cmpd_ids,
                                  ref_subnodes,
                                  ref_synthons,
                                  args.output_dir,
                                  args.working_dir,
                                  args.scoring_dir,
                                  args.type_merge,
                                  args.substructure_dir,
                                  similarities,
                                  used_subnodes,
                                  used_synthons,
                                  original_similarities,
                                  original_names,
                                  rev_synthons,
                                  intermed_nodes1,
                                  intermed_nodes2,
                                  limit,
                                  minimize_limit,
                                  args.n_cpus,
                                  args.move_files)

    pipeline.run_filtering()
    pipeline.run_scoring()


if __name__ == "__main__":
    main()
