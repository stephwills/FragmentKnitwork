# Fragment Knitwork

![Alt text](imgs/knitwork.png)

This repository contains the code for running fragment merging using the new fragment network pipeline.
The pipeline 'knits' together substructures from crystallographic fragment hits, creating merge and linker-like follow-up compounds.
This is an updated version of the fragment_network_merges repository and addresses some of the limitations of the original pipeline and provides some additional feature enhancements
(for example, querying is done at the substructure level and the database can search for pharmacophorically similar substructures to incorporate).
The method uses a neo4j graph database, constructed in the same way as described in the original Astex implementation (Hall et al. *J. Med. Chem.* **2017** *60*(14), 6440-6450).

## Overview
There are three main modules:
1) **Fragment:** This module enumerates possible pairs substructures from the crystallographic fragments, and performs some filtering of those pairs to identify the best opportunities.
2) **Knitwork**: This module 'knits' together the fragments, i.e. takes care of all the database querying to the neo4j database and retrieves possible merged molecules.
3) **Quilter**: this if the f(qu)iltering part of the pipeline and provides functions for generating conformers, largely using RDKit's ph4 implementation and Fragmenstein (Ferla et al. *ChemRxiv* **2024** DOI:10.26434/chemrxiv-2024-17w01). There are options for using PyRosetta for minimization or not using PyRosetta (as it requires an academic license).

The utils folder contains all the utility functions for each module.
There are several parameters throughout the pipeline that can be controlled either through command line arguments or can be specified in the Config files in the utils directory. 

## Installation
The following packages are used in this repo:

- RDKit
- fragmenstein (https://github.com/matteoferla/Fragmenstein)
- tqdm
- joblib
- pymol-open-source
- MDAnalysis
- ProLIF
- neo4j
- PyRosetta (requires academic license)

An example of how to do this is shown below. For PyRosetta installation, replace the username and password as appropriate.

```
conda create --name fragnet
conda activate fragnet -y
conda install -c conda-forge rdkit -y
conda install -c https://username:password@west.rosettacommons.org/pyrosetta/conda/release/ pyrosetta -y
pip install tqdm joblib fragmenstein MDanalysis prolif
git clone https://github.com/schrodinger/pymol-open-source.git && cd pymol-open-source && python setup.py install && cd .. && rm -r pymol-open-source
```

## Docker container for running the neo4j database

https://github.com/stephwills/docker-neo4j contains the instructions for how to run the neo4j database using a Docker container.
The database files are not currently available but will be added shortly (plus a small test version).

## Usage

### Fragment

This module enumerates compatible pairs of substructures for merging, filtering them according to overlap and checking whether they make interactions with the protein. The initial enumeration involves querying the dataabase to retrieve constituent substructures for a fragment node.

To run the enumeration, use the file `FragmentKnitwork/Fragment/runEnumeration.py`. This can be run on the command line.
You can provide one of two types of inputs: either a comma-delimited text file with each line proving a name for the fragment, the fragment SMILES and the path to a .mol file for the fragment.
You must be connected to the neo4j database for this module to work.

The `--record_equiv_synthon` and `--r_group_expansions` flags are used to generate files that are used later in the bioisosteric merging queries and for recapitulating R groups seen in the fragments.

```angular2html
usage: runEnumeration.py [-h] [--input_txt_file INPUT_TXT_FILE]
                         [--fragment_names FRAGMENT_NAMES [FRAGMENT_NAMES ...]]
                         [--target TARGET] [--fragalysis_formatted]
                         [--write_smiles_to_file] [--record_equiv_synthon]
                         [--r_group_expansions]
                         [--pairs_to_ignore_json PAIRS_TO_IGNORE_JSON]
                         [--output_dir OUTPUT_DIR]

options:
  -h, --help            show this help message and exit
  --input_txt_file INPUT_TXT_FILE
                        comma delimited txt file with names, SMILES columns, mol files
                        and apo-desolv pdb files
  --fragment_names FRAGMENT_NAMES [FRAGMENT_NAMES ...]
                        provide if data is Fragalysis-formatted
  --target TARGET       target name for Fragalysis-formatted data
  --fragalysis_formatted
                        flag to indicate that data is Fragalysis-formatted
  --write_smiles_to_file
                        write out SMILES pairs to a file
  --record_equiv_synthon
                        (recommended) record equivalent synthons/subnodes for fragments,
                        used for some of the more complex queries
  --r_group_expansions  (recommended) record possible r group expansions for future
                        queries
  --pairs_to_ignore_json PAIRS_TO_IGNORE_JSON
                        pairs to not run enumeration before, provide in format (remember
                        pairs are asymmetric) [["x0000_2A", "x000_1A"], ["x0000_1A",
                        "x000_2A"]
  --output_dir OUTPUT_DIR
                        where to save the output data
```

**Example usage when using an input txt file:**
```
python Fragment/runEnumeration.py --input_txt_file fnames.txt --target TARGET
--record_equiv_synthon --r_group_expansions --output_dir enumeration_output/
```

Alternatively, if the files are formatted in the old Fragalysis download format (i.e. named like x0072_0A), you can simply provide a list of names of the fragments you want to run the pipeline on.
See `FragmentKnitwork/Fragment/README.md` for more details.
Note: The following steps currently only support data in the Fragalysis format described! 

**Example usage for fragalysis formatted data:**
```
python Fragment/runEnumeration.py --fragment_names x0234_0A x0123_0A x0756_0A --target TARGET --fragalysis_formatted --record_equiv_synthon --r_group_expansions --output_dir enumeration_output/
```

There are several parameters controlling how strict the enumeration is, these can be set in `FragmentKnitwork/utils/fragmentConfig.py`. 

### Knitwork

This module allows us to run the queries within the neo4j database. There are two main types of query: for identifying perfect merges (whereby exact substructures are merged into a single compound) or bioisosteric merges (whereby a substructure is replaced with one that is 'bioisosterically' similar through the use of a similarity calculation between pharmacophore fingerprints). 

The queries can be run using the runKnitting function defined in `Knitwork/runKnitting.py` and can also be run via the command line.
The code takes as input the `substructure_pair.json` generated by `FragmentKnitwork/Fragment/runEnumeration.py`. 

```angular2html
usage: runKnitting.py [-h] [--substructure_pair_file SUBSTRUCTURE_PAIR_FILE] [--descriptor {prop_pharmfp}] [--pure_search]
                      [--working_dir WORKING_DIR] [--output_dir OUTPUT_DIR] [--n_parallel N_PARALLEL] [--run_parallel]
                      [--clear_working_dir] [--limit LIMIT] [--prolif_prioritization] [--target TARGET]
                      [--substructure_dir SUBSTRUCTURE_DIR] [--max_prioritize MAX_PRIORITIZE] [--r_group_search R_GROUP_SEARCH]

options:
  -h, --help            show this help message and exit
  --substructure_pair_file SUBSTRUCTURE_PAIR_FILE
                        substructure_pairs.json file generated by runEnumeration.py, if not provided, checks within substructure_dir
  --descriptor {prop_pharmfp}
                        descriptor property for bioisosteric merging search, currently only supports prop_pharmfp!
  --pure_search         flag for if you want to search for only pure merges
  --working_dir WORKING_DIR
                        where to save intermediate files
  --output_dir OUTPUT_DIR
                        where to save output files
  --n_parallel N_PARALLEL
                        number of parallel queries to run
  --run_parallel        whether to run queries in parallel
  --clear_working_dir   whether to clear the contents of the working dir
  --limit LIMIT         limit to the number of compounds to retrieve per query
  --prolif_prioritization
                        whether to prioritize the substructures used for those that make an interactions
  --target TARGET       name of the target (in fragalysis format)
  --substructure_dir SUBSTRUCTURE_DIR
                        location of the enumeration dir creater by runEnumeration.py
  --max_prioritize MAX_PRIORITIZE
                        if prioritizing, the max number of compounds to prioritize
  --r_group_search R_GROUP_SEARCH
                        whether you want to run an r group search on the compounds automatically using R groups from the original
                        fragments
```
**Example usage for perfect merging:**

- To run a search for perfect merges, use the `--pure_search` flag. Queries can be run in parallel by specifying `--run_parallel` and providing the number of queries to run in parallel with `--n_parallel`.
A limit can be imposed on the number of compounds to be found from the search for each substructure pair using `--limit`.

```
python Knitwork/runKnitting.py --substructure_pair_file enumeration_dir/substructure_pairs.json --pure_search --working_dir work_dir --output_dir out_dir --run_parallel --n_parallel 2 --limit 100 --clear_working_dir
```

**Example usage for bioisosteric merging:**

- To run a bioisosteric search (this is a different query so won't overlap with the perfect merge results), run the command without the `--pure_search` flag. The descriptor to be used to find replacement substructures is specified using `--descriptor` So far, only the pharmacophore fingerprint descriptor (`prop_pharmfp`) is supported for finding replacement substructures.

```
python Knitwork/runKnitting.py --substructure_pair_file enumeration_dir/substructure_pairs.json --descriptor prop_pharmfp --working_dir work_dir --output_dir out_dir --run_parallel --n_parallel 2 --limit 100 --clear_working_dir
```

**Example usage for bioisosteric merging with prioritization:**

- We can optionally prioritize the results of the bioisosteric search by only keeping replacement substructures that potentially recapitulate an interaction made by the original substructure. We can also order them by number of interactions made and take the top N. Currently, this is only supported for data that is in the Fragalysis format and will require you to specify the substructure dir created in `Fragment/runEnumeration.py`, as it contains files containing the substructures in their original conformations. This is done as follows:

```
python Knitwork/runKnitting.py --substructure_pair_file enumeration_dir/substructure_pairs.json --descriptor prop_pharmfp --working_dir work_dir --output_dir out_dir --run_parallel --n_parallel 2 --limit 100 --clear_working_dir --prolif_prioritization --max_prioritize 50 --target TARGET --substructure_dir subs_dir
```

- We can also run an R-group expansion on the perfect and bioisosteric merges using the `--r_group_search` flag (this relies on you using the revelant flags during the runEnumeration stage). This is largely untested (for the paper, the R group expansion was done after alignment of the placed molecules using the RGroupExpansionOnScored.py file, to check whether the scaffolds can be placed first).

**Example usage for running a second round of bioisosteric merging:**

- We can run an additional bioisosteric replacement for merges (replacing the seed node substructure in the query) using `FragmentKnitwork/Knitwork/reverseQueries.py`. This uses the scored sdf file for the bioisosteric merges generated by Quilter (see below).
- We can specify the max number of bioisosteric replacements for each placed merge with `--limit_results`
- The `--linker_option` refers to whether we want to allow no change in the linker (strict), a change in the linker (loose), or to run both options (both)\

```
python Knitwork/reverseQueries.py --sdf_file scored_dir/scored_mols.sdf --substructure_dir enumeration_dir/ --working_dir work_dir --output_dir out_dir --n_parallel 2 --limit_results 20 --threshold 0.9 --linker_option strict
```

### Quilter

Qu(f)ilter is used to filter the merges by generating conformations and using a few simple filters to remove unfavourable compounds. There are two options: to run filtering with minimization on ALL merges (using Fragmenstein's Victor) or by running filtering without minimization on all compounds and then performing filtering again with minimization on the top X. This is done by running the `Quilter/runQuilterPipeline.py` file.

*Code has been recently refactored to make it easier to use -- more thorough tests will be added in test/ to ensure all is working OK.*

```angular2html
usage: runQuilterPipeline.py [-h] --json_file JSON_FILE --pdb_file PDB_FILE [--fragmentA FRAGMENTA] [--fragmentB FRAGMENTB]
                             [--target TARGET] [--substructure_dir SUBSTRUCTURE_DIR] [--output_dir OUTPUT_DIR] [--working_dir WORKING_DIR]
                             [--scoring_dir SCORING_DIR] [--n_cpus N_CPUS] [--type_merge {pure_merge,impure_merge,reverse_merge}]
                             [--limit_num_run] [--limit_num_minimize] [--max_num_filter MAX_NUM_FILTER]
                             [--max_num_minimize MAX_NUM_MINIMIZE] [--move_files] [--use_wictor]

options:
  -h, --help            show this help message and exit
  --json_file JSON_FILE
                        This should be one of the json files generated for a fragment pair by runKnitting.py
  --pdb_file PDB_FILE   apo-desolv pdb file to use for placing the merges
  --fragmentA FRAGMENTA
                        The name of the first fragment in the pair (how it is named in the Fragalysis format)
  --fragmentB FRAGMENTB
                        The name of the second fragment in the pair (how it is named in the Fragalysis format)
  --target TARGET       The name of the target (how it is named in the Fragalysis format)
  --substructure_dir SUBSTRUCTURE_DIR
                        The dir generated by runEnumeration.py
  --output_dir OUTPUT_DIR
                        where to save output data
  --working_dir WORKING_DIR
                        where to save intermediate files
  --scoring_dir SCORING_DIR
                        where to save the scored ligands
  --n_cpus N_CPUS       the number of CPUs to parallelise over
  --type_merge {pure_merge,impure_merge,reverse_merge}
                        The type of merge, a pure (perfect) merge, or a bioisosteric merge where one (impure merge) or two substructures
                        (reverse merge) has been replaced
  --limit_num_run       whether to limit the number of compounds overall filtered
  --limit_num_minimize  only relevant if using wictor pre-filter, whether to limit the number of compounds minimized (after the first pre-
                        filter)
  --max_num_filter MAX_NUM_FILTER
                        the max number of compounds to filter overall
  --max_num_minimize MAX_NUM_MINIMIZE
                        only relevant if using wictor pre-filter, the max number of compounds to minimize
  --move_files          (recommended) whether to move files back to the output directory
  --use_wictor          whether to use the wictor pre-filter (without minimization) then prioritize the top X
```
**Example usage without pre-filter**

Here the `--limit_num_run` flag indicates we want to limit the number of molecules run through the pipeline, which is controlled by `--max_num_filter`. 

```
python Quilter/runQuilterPipeline.py --json_file x0000_1A-x0000_2A_prop_pharmfp_impure_merges.json --pdb_file x0000_1A_apo-desolv.pdb --fragmentA x0000_1A --fragmentB x0000_2A --target TARGET --substructure_dir enumeration_dir --output_dir out_dir --working_dir work_dir --scoring_dir scored_dir --n_cpus 4 --type_merge impure_merge --limit_num_run --max_num_filter 1000 --move_files
```

**Example usage with pre-filter**

Here we can specify that we want to use the pre-filter pipeline with Wictor with `--use_wictor`. We control the number of compounds we want to be subsequently minimized with the flag `--limit_num_minimize` and setting the number with `--max_num_minimize`.

```
python Quilter/runQuilterPipeline.py --json_file x0000_1A-x0000_2A_prop_pharmfp_impure_merges.json --pdb_file x0000_1A_apo-desolv.pdb --fragmentA x0000_1A --fragmentB x0000_2A --target TARGET --substructure_dir enumeration_dir --output_dir out_dir --working_dir work_dir --scoring_dir scored_dir --n_cpus 4 --type_merge impure_merge --limit_num_run --limit_num_minimize --max_num_filter 1000 --max_num_minimize 500 --move_files --use_wictor
```

Note: in the paper, R-group expansion was run on the filtered compounds to limit computational time (rather than at the start). The code aligning these is in `runRGroupExpansionAlignment.py`.

### Tests

There are some limited tests available in `test/`. 