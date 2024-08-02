import os
import subprocess
from argparse import ArgumentParser

from FragmentKnitwork.utils import quilterConfig as config
from FragmentKnitwork.utils.utils import get_protein


def main():
    parser = ArgumentParser()
    parser.add_argument('--py_file_path', help='path to runAlignment.py file')
    parser.add_argument('--json_dir')
    parser.add_argument('--target', required=False)
    parser.add_argument('--method_tag', help='tag used in naming the json files')
    parser.add_argument('--type_merge', choices=['pure_merge', 'impure_merge'])
    parser.add_argument('--substructure_dir', required=False, default=config.SUBSTRUCTURE_DIR)
    parser.add_argument('--output_dir', required=False, default=config.OUTPUT_DIR)
    parser.add_argument('--working_dir', required=False, default=config.WORKING_DIR)
    parser.add_argument('--n_cpus', type=int, required=False, default=config.N_CPUS)
    parser.add_argument('--parallel', action='store_true', default=config.PARALLEL)
    parser.add_argument('--min_files', action='store_true', default=config.MIN_FILES)
    parser.add_argument('--move_files', action='store_true', default=config.MOVE_FILES)
    parser.add_argument('--limit_num_run', action='store_true')
    parser.add_argument('--max_num_filter', type=int, required=False, default=None)
    parser.add_argument('--max_num_minimize', type=int, required=False, default=None)
    parser.add_argument('--run_scoring', action='store_true')
    args = parser.parse_args()

    files = [file for file in os.listdir(args.json_dir)]
    pairs = [file.replace(f"_{args.method_tag}", '').replace('.json', '') for file in files]
    files = [os.path.join(args.json_dir, file) for file in files]
    fAs = [pair.split('-')[0] for pair in pairs]
    fBs = [pair.split('-')[1] for pair in pairs]

    pdb_files = [get_protein(args.target, fA) for fA in fAs]

    print(len(files), 'files ready for alignment')

    for file, fA, fB, pair, pdb_file in zip(files, fAs, fBs, pairs, pdb_files):
        print('Running alignment for pair:', pair)
        command = f"""python {args.py_file_path} --json_file {file} --pdb_file {pdb_file} --target {args.target} --fragmentA {fA} --fragmentB {fB} --substructure_dir {args.substructure_dir} --output_dir {args.output_dir} --working_dir {args.working_dir} --n_cpus {args.n_cpus} --type_merge {args.type_merge} --parallel --min_files --move_files --limit_num_run --max_num_filter {args.max_num_filter} --max_num_minimize {args.max_num_minimize}"""
        if args.run_scoring:
            command = command + ' --run_scoring'
        print('Command:')
        print(command)
        subprocess.run(command, shell=True)


if __name__ == "__main__":
    main()