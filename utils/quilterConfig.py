########## config for running alignment with fragmenstein ##########
ALIGNMENT_FDEF = 'AlignmentFeatures.fdef'                       # these are the features used for aligning the merge to the original fragment substructure ph4s
PH4_RADII = 0.5                                                 # default ph4 radii
ALIGNMENT_SCORE = 'SuCOS'                                       # how to score the alignment (only supports highest SuCOS right now)
SCORING_MODE = 'max'                                            # how to select for best score (i.e. max score when evaluating sucos)
EMBEDDINGS_COUNT = 10                                           # number of embeddings to generate in ph4 embedding
N_FEAT_COMBINATIONS = 100                                       # number of possible pharmacophore feature combinations to try
ATOM_MATCH_THRESHOLD = 1.0                                      # distance to use for selecting atoms to map to each other when generating custom map
PYROSETTA_MINIMIZATION = True                                   # whether to perform minimization using PyRosetta
DDG_FILTER = True                                               # whether to apply DDG filter for fragmenstein
RMSD_THRESHOLD = 2.0                                            # RMSD threshold for placement against parent substructures
REMOVE_WICTOR = True                                            # whether to remove intermediate Wictor directories
COVALENT_RESI = '30A'                                           # covalent resi to set in Fragmenstein


########## config for overall running alignment and scoring ##########
OUTPUT_DIR = '/home/swills/test/output'                         # output dir for alignment files
WORKING_DIR = '/home/swills/test/working'                       # working dir for alignment files
SUBSTRUCTURE_DIR = '/home/swills/test_substructures'            # where substructures were saved (generated in Fragment/run_enumeration)
PARALLEL = True                                                 # whether to run in parallel
N_CPUS = 2                                                      # number of CPUs to parallelise with
MIN_FILES = True                                                # whether to only keep the minimal alignment files needed (rather than all Fragmenstein output)
MOVE_FILES = True                                               # whether to move files to output directory at the end (from working)
INFO_FILE_EXISTS = True                                         # whether the substructure_pair file has been generated
LIMIT_NUM_RUN = False                                           # whether to limit the number of alignment molecules
MAX_NUM_FILTER = 500                                            # if LIMIT_NUM_RUN, how many to limit to


########## config for specific filters ##########
RUN_OVERLAP_CHECK = True                                        # whether to check overlap between mol and substructures
RUN_ENERGY_CHECK = True                                         # whether to check energy of the generated conformation
ENERGY_CHECK_N_CONF = 50                                        # how many conformers to generate when running the energy check
ENERGY_CHECK_THRESHOLD = 7.0                                    # threshold for energy check (ratio of avg unconstrained against constrained conformation)


########## config for interaction calculations ##########
INTERACTIONS = [
    "hydrophobic_contacts",
    "pi_stacking",
    "pi_cation",
    "hbond_pdon",
    "hbond_ldon",
    "saltbridge_lneg",
    "saltbridge_pneg",
    "saltbridge_metal",
    "halogen",
]
PROLIF_INTERACTIONS = [
    'Hydrophobic',
    'HBAcceptor',
    'HBDonor',
    'Cationic',
    'Anionic',
    'CationPi',
    'PiCation',
    'PiStacking'
]