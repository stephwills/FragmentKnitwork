import os
import shutil
from argparse import ArgumentParser

from FragmentKnitwork.utils.quilterUtils import calc_prolif_interaction
from FragmentKnitwork.utils.utils import get_mol, get_protein, dump_json


def get_all_fragment_fps(fragments, working_dir, output_dir, target):
    """
    Get fragment IFPs with ProLIF for use in scoring later on

    :param fragments:
    :param working_dir:
    :param output_dir:
    :param target:
    :return:
    """
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)

    output_fname = os.path.join(output_dir, 'fragment_interactions.json')
    if os.path.exists(output_fname):
        print(output_fname, 'already exists')
        return None

    interaction_res = {}

    for frag in fragments:
        mol_file = get_mol(target, frag)
        apo_file = get_protein(target, frag, protonated=True, desolv=True)
        ints = calc_prolif_interaction(mol_file, apo_file, lig_protonated=False, write_interactions_file=False)
        interaction_res[frag] = ints

    dump_json(interaction_res, output_fname)
    shutil.rmtree(working_dir)


def main():
    parser = ArgumentParser()
    parser.add_argument('--fragment_names', nargs="+")
    parser.add_argument('--target')
    parser.add_argument('--output_dir')
    parser.add_argument('--working_dir')
    args = parser.parse_args()

    if not os.path.exists(args.working_dir):
        os.mkdir(args.working_dir)

    output_fname = os.path.join(args.output_dir, 'fragment_interactions.json')
    if os.path.exists(output_fname):
        print(output_fname, 'already exists')
        return None

    interaction_res = {}

    fragments = args.fragment_names
    for frag in fragments:
        mol_file = get_mol(args.target, frag)
        apo_file = get_protein(args.target, frag, protonated=True, desolv=True)
        ints = calc_prolif_interaction(mol_file, apo_file, lig_protonated=False, write_interactions_file=False)
        interaction_res[frag] = ints

    dump_json(interaction_res, output_fname)


if __name__ == "__main__":
    main()
