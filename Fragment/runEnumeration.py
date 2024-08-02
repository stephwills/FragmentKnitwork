from argparse import ArgumentParser
import os

from FragmentKnitwork.Fragment.enumerate import enumeration


def runEnumeration(fragment_names, fragment_smiles, target, mol_files, pdb_files, output_dir, ignore_pairs=None,
                   write_smiles_to_file=True, r_group_expansions=True, record_equiv_synthon=True):
    """

    :param fragment_names:
    :param fragment_smiles:
    :param target:
    :param mol_files:
    :param pdb_files:
    :param output_dir:
    :param ignore_pairs:
    :param write_smiles_to_file:
    :param r_group_expansions:
    :param record_equiv_synthon:
    :return:
    """
    # run enumeration for the substructures
    substructure_pair_fname = enumeration(fragment_names, fragment_smiles, target, mol_files=mol_files, pdb_files=pdb_files, substructure_dir=output_dir,
                                          ignore_pairs=ignore_pairs)
    print('Enumeration finished')

    if write_smiles_to_file:
        # write list of smiles pairs to file -- useful if we want to get an equivalent list for similarity search later
        from FragmentKnitwork.utils.utils import get_smiles, load_json
        substructure_pair_data = load_json(substructure_pair_fname)
        fragment_pairs = set()
        for fragment_pair in substructure_pair_data:
            p = fragment_pair.split('-')
            p.sort()
            fragment_pairs.add(f"{get_smiles(target, p[0])},{get_smiles(target, p[1])}")
        fragment_pairs = list(fragment_pairs)
        fragment_pairs.sort()
        txt_file = os.path.join(output_dir, f'{target}_smiles_pairs.txt')
        with open(txt_file, 'w') as f:
            for line in fragment_pairs:
                f.write(line)
                f.write('\n')

    if r_group_expansions:
        # this retrieves the possible r groups that can be used for r group expansion queries later on
        r_group_expansion_data = {}
        from FragmentKnitwork.utils.fragmentUtils import run_r_group_expansion
        from FragmentKnitwork.utils.utils import dump_json
        for fragment, smi in zip(fragment_names, fragment_smiles):
            r_group_expansions = run_r_group_expansion(smi)
            r_group_expansion_data[fragment] = r_group_expansions
        dump_json(r_group_expansion_data, os.path.join(output_dir, 'r_group_expansions.json'))

    if record_equiv_synthon:
        # this is useful to record for some of the more complex queries later on (if we want to replace the seed node in the query)
        from FragmentKnitwork.Fragment.substructure import equivalent_synthon_subnodes
        equiv_json = os.path.join(output_dir, 'equivalent_subnodes.json')
        equivalent_synthon_subnodes(fragment_names, fragment_smiles, saveJson=equiv_json)


def main():
    parser = ArgumentParser()
    parser.add_argument('--input_txt_file', required=False, help='comma delimited txt file with names, SMILES columns, mol files and apo-desolv pdb files')
    parser.add_argument('--fragment_names', nargs="+", required=False, help='provide if data is Fragalysis-formatted')
    parser.add_argument('--target', required=False, help='target name for Fragalysis-formatted data')
    parser.add_argument('--fragalysis_formatted', action='store_true', help='flag to indicate that data is Fragalysis-formatted')
    parser.add_argument('--write_smiles_to_file', action='store_true', help='write out SMILES pairs to a file ')
    parser.add_argument('--record_equiv_synthon', action='store_true', help='(recommended) record equivalent synthons/subnodes for fragments, used for some of the more complex queries')
    parser.add_argument('--r_group_expansions', action='store_true', help='(recommended) record possible r group expansions for future queries')
    parser.add_argument('--pairs_to_ignore_json', default=None, help='pairs to not run enumeration before, provide in format (remember pairs are asymmetric) [["x0000_2A", "x000_1A"], ["x0000_1A", "x000_2A"]')
    parser.add_argument('--output_dir', help="where to save the output data")
    args = parser.parse_args()

    # read data provided -> get the fragment names, smiles and mol files
    if args.input_txt_file:
        import csv
        fragment_names = []
        fragment_smiles = []
        mol_files = []
        pdb_files = []
        with open(args.input_txt_file, 'r') as f:
            for row in csv.reader(f, delimiter=','):
                fragment_names.append(row[0])
                fragment_smiles.append(row[1])
                mol_files.append(row[2])
                pdb_files.append(row[3])

    if args.fragalysis_formatted:
        from FragmentKnitwork.utils.utils import get_mol, get_smiles, get_protein
        if not args.fragment_names or not args.target:
            raise ValueError('Please provide fragment names and target if data is fraglaysis formatted')
        fragment_names = args.fragment_names
        mol_files = [get_mol(args.target, name) for name in fragment_names]
        pdb_files = [get_protein(args.target, name, protonated=True, desolv=True) for name in fragment_names]
        fragment_smiles = [get_smiles(args.target, name) for name in fragment_names]

    # check if there are any pairs to be ignored (e.g. if you want to focus on pairs between subsites)
    print(len(fragment_names), 'fragments provided for enumeration')
    if args.pairs_to_ignore_json:
        from FragmentKnitwork.utils.utils import load_json
        pairs_to_ignore = load_json(args.pairs_to_ignore_json)
    else:
        pairs_to_ignore=None

    runEnumeration(fragment_names, fragment_smiles, args.target, mol_files, pdb_files, args.output_dir, ignore_pairs=pairs_to_ignore,
                   write_smiles_to_file=args.write_smiles_to_file, r_group_expansions=args.r_group_expansions, record_equiv_synthon=args.record_equiv_synthon)


if __name__ == "__main__":
    main()
