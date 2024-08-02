# this is a customized fdef file used to calculate pharmacophore features for substructures
# in the paper version, aliphatic rings have been added (to preserved shape)
# and xenon atoms (which denote attachment points)
FINGERPRINT_FDEF = 'FeatureswAliphaticXenon.fdef'
FINGERPRINT_MAXPOINTCOUNT = 2                                       # max point count for rdkit pharm fp calculation
FINGERPRINT_BINS = [(0, 2), (2, 5), (5, 8)]                         # fingerprint binds for rdkit pharm fp calculation

SINGLE_EXPANSION_LIMIT = None                                       # limit to the number of R-group query results we want to retrieve
RESULTS_LIMIT = None                                                # limit the number of results for the pure/impure queries

NUM_HOPS = 2                                                        # number of optional hops (expansions) to be made before the substructure expansion
INCLUDE_PURE = False
SIMILARITY_THRESHOLD = 0.9                                          # threshold for similarity for selecting a replacement substructure
SIMILARITY_METRIC = 'usersimilarity.tanimoto_similarity'            # name of the similarity metric to use within neo4j
DESCRIPTOR_NAME = 'prop_pharmfp'                                    # the name added to the edges for the descriptors describing substructures

REVERSE_QUERY_LIMIT = 5                                             # limit the number of results from the reverse query
REVERSE_QUERY_STRICT = True                                         # whether the reverse query should be strict (no changes made to the linker) or not

OUTPUT_DIR = None
WORKING_DIR = None