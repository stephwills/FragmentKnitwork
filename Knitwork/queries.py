from FragmentKnitwork.utils import knitworkConfig as config
from FragmentKnitwork.utils.dbUtils import driver


def single_expansion(smiles: str, query_synthon: str, limit=config.SINGLE_EXPANSION_LIMIT):
    """
    Simple query used for R-group expansion if we have an exact substructure we want to expand with. There
    are no intermediate hops -> only direct expansions off the molecule.

    :param smiles: smiles of the molecule we want to expand
    :param query_synthon: smiles for the substructure we want to add, e.g. '[Xe]Cl'
    :param limit: optional to limit the number of results we want to retrieve
    :return: a list of the expansions SMILES and associated compound IDs
    """
    query = (
        """
        MATCH (a:F2 {smiles: $smiles})<-[e:FRAG]-(c:Mol)
        WHERE e.prop_synthon=$query_synthon
        WITH c.smiles as smi, c.cmpd_ids as ids
        RETURN smi, ids
        """
        % {"smiles": smiles,
           "query_synthon": query_synthon}
    )
    if limit:
        query = query + f" LIMIT {limit}"

    expansions = []
    cmpd_ids = []

    with driver.session() as session:
        results = session.run(query, smiles=smiles, query_synthon=query_synthon)
        for res in results:
            if res['smi'] not in expansions:
                expansions.append(res['smi'])
                cmpd_ids.append(res['ids'])
    return expansions, cmpd_ids


def pure_expansion(smiles: str, query_synthon: str, num_hops: int = config.NUM_HOPS, limit=config.RESULTS_LIMIT):
    """
    Query for running a pure expansion (whether the exact expansion substructure is specified, i.e. no analogues
    are incorporated). Starts from a specified node smiles, a series of optional hops (which are specifically
    expansions) and a final expansion, where the specified substructure is added.

    :param smiles: SMILES of the node to start from
    :param query_synthon: SMILES of the substructure to incorporate (e.g. '[Xe]c1ccccc1')
    :param num_hops: number of optional hops (i.e. up to two expansions can be made)
    :param limit: number of results to limit the query to (e.g. only retrieve up to 500 possible expansions)
    :return: list of SMILES of the exapnsions and associated compound IDs
    """
    query = (
        """
        MATCH (a:F2 {smiles: $smiles})<-[:FRAG*0..%(num_hops)d]-(b:F2)<-[e:FRAG]-(c:Mol)
        WHERE e.prop_synthon=$query_synthon
        WITH c.smiles as smi, c.cmpd_ids as ids
        RETURN smi, ids
        """
        % {"smiles": smiles,
           "query_synthon": query_synthon,
           "num_hops": num_hops}
    )
    if limit:
        query = query + f" LIMIT {limit}"
    expansions = []
    cmpd_ids = []
    with driver.session() as session:
        results = session.run(query, smiles=smiles, query_synthon=query_synthon, num_hops=num_hops)
        for res in results:
            expansions.append(res['smi'])
            cmpd_ids.append(res['ids'])

    return expansions, cmpd_ids


def impure_expansion(smiles: str, vector, query_synthon: str, metric: str = config.SIMILARITY_METRIC, desc: str = config.DESCRIPTOR_NAME,
                     threshold: float = config.SIMILARITY_THRESHOLD, num_hops=config.NUM_HOPS, limit=config.RESULTS_LIMIT):
    """
    Query for identifying the bioisosteric/impure merges. Starts from a specified node SMILES, makes up to a specified
    number of optional hops (expansions), and a final expansion where a substructure is incorporated where the substructure
    has similarity above a certain threshold to a query substructure (by supplying the vector for the query substructure).
    Similarity is calculated on the fly using a neo4j function (may have to be manually created and added to the database).

    :param smiles: SMILES of the initial node
    :param vector: pre-calculated fp (or whatever descriptor) for the query substructure for calculating similarity against
    :param query_synthon: SMILES of the query substructure (so we don't retrieve the exact structure)
    :param metric: NOT USED RIGHT NOW, name of the similarity function within the database, couldn't get this to work so has been manually specified
    :param desc: name of the descriptor (just prop_pharmfp right now) to specify which metric/query type to use
    :param threshold: similarity threshold for selecting substructure analogues
    :param num_hops: number of optional hops to make before expansion
    :param limit: an optional limit the limit the number of results
    :return: lists of the expansion SMILES, actual substructures incorporated, similarity and associated compound IDs
    """
    # here we specifically try to find merges that don't incorporate the exact substructure (can adapt this later if needed)
    if desc == 'prop_pharmfp':
        query = (
            """
            MATCH (a:F2 {smiles: $smiles})<-[:FRAG*0..%(num_hops)d]-(b:F2)<-[e:FRAG]-(c:Mol)
            WHERE EXISTS(e.prop_pharmfp)
            WITH usersimilarity.tanimoto_similarity(e.prop_pharmfp, $vector) as sim, c.smiles as smi, e.prop_synthon as syn, c.cmpd_ids as ids
            WHERE sim >= $threshold
            AND NOT e.prop_synthon=$query_synthon
            RETURN smi, syn, sim, ids
            """
            % {"smiles": smiles,
               "vector": vector,  # a pre-computed fp for the query substructure used to calculate similarity
               "threshold": threshold,  # a threshold for similarity to identify replacement subtructures
               "query_synthon": query_synthon,  # the SMILES for the query synthon (to make sure we don't retrieve exact expansions with it)
               "metric": metric,  # TODO: not used here: we have manually specified the name of the similarity function
               "num_hops": num_hops}  # the max number of optional hops (expansions) to make before incorporating the substructure
        )
    # TODO: this is not used in the paper as USRCAT requires 3D conformer generation
    if desc == 'prop_usrcat':
        query = (
            """
            MATCH (a:F2 {smiles: $smiles})<-[:FRAG*0..%(num_hops)d]-(b:F2)<-[e:FRAG]-(c:Mol)
            WHERE EXISTS(e.prop_usrcat)
            WITH usersimilarity.usrcat_similarity_xenon(e.prop_usrcat, $vector) as sim, c.smiles as smi, e.prop_synthon as syn, c.cmpd_ids as ids
            WHERE sim >= $threshold
            AND NOT e.prop_synthon=$query_synthon
            RETURN smi, syn, sim, ids
            """
            % {"smiles": smiles,
               "vector": vector,
               "threshold": threshold,
               "query_synthon": query_synthon,
               "metric": metric,
               "num_hops": num_hops}
        )
    if limit:
        query = query + f" LIMIT {limit}"
    expansions = []
    synthons = []
    similarities = []
    cmpd_ids = []
    with driver.session() as session:
        results = session.run(query, smiles=smiles, vector=vector, threshold=threshold,
                              query_synthon=query_synthon, metric=metric, num_hops=num_hops)
        for res in results:
            expansions.append(res['smi'])
            synthons.append(res['syn'])
            similarities.append(res['sim'])
            cmpd_ids.append(res['ids'])

    return expansions, synthons, similarities, cmpd_ids


def replace_substructure(smiles: str, synthon: str, threshold: float = config.SIMILARITY_THRESHOLD,
                         rev_limit=config.REVERSE_QUERY_LIMIT, strict: bool = config.REVERSE_QUERY_STRICT,
                         first_synthon=None):
    """
    This query is used when we have a merge (after running the impure_expansion query) and want to go back and replace
    the original substructure represented by the seed node in the query, termed a 'reverse query'. Can run a strict
    query (where the linker doesn't change at all) or a non-strict query (where it can undergo an un-specified contraction
    and expansion).

    :param smiles: SMILES of the initial node to start from
    :param synthon: synthon of the substructure we want to replace (representing the initial node in first query)
    :param threshold: threshold for similarity for substructure
    :param rev_limit: optional limit to limit the number of results
    :param strict: bool of whether to run strict query or not
    :param first_synthon: the substructure that was added in the first query (that we want to preserve)
    :return: lists of expansions, replacement substructures added, similarty values, compound ids, and
    intermediate nodes (one or two lists, second empty depending on whether strict or not)

    """
    # strict queries are if we want to make NO CHANGE to the linker (so we only remove the substructure, and add a replacement)
    if strict:
        query = (
            """
            MATCH (a:F2 {smiles: $smiles})-[e1:FRAG]->(b:F2)<-[e2:FRAG]-(c:Mol)
            WHERE e1.prop_synthon=$synthon
            AND EXISTS(e2.prop_pharmfp)
            WITH usersimilarity.tanimoto_similarity(e1.prop_pharmfp,e2.prop_pharmfp) as sim, e2.prop_synthon as syn, c.smiles as smi, c.cmpd_ids as ids, b.smiles as intermed1
            WHERE sim >= $threshold
            AND NOT e2.prop_synthon=$synthon
            RETURN sim, syn, smi, ids, intermed1 LIMIT $rev_limit
            """
            % {"smiles": smiles,
               "synthon": synthon,
               "threshold": threshold,
               "rev_limit": rev_limit}
        )
    # non-strict query where we can make a change in the linker (an unspecified contraction followed by unspecified
    # expansion. Several controls have to be specified: 1) we don't want to remove and add back the same thing in the linker,
    # 2) we don't want to remove the OTHER substructure initially added in the first query (impure_merge), and 3)
    # we don't want to add back the substructure we're replacing in the final expansion (we're looking for analogs)
    else:
        query = (
            # TODO: How to specify with 'up to one hop' the synthon
            # AND e2.prop_synthon=$first_synthon
            """ 
            MATCH (a:F2 {smiles: $smiles})-[e1:FRAG]->(b:F2)-[e2:FRAG]->(c:F2)<-[e3:FRAG]-(d:F2)<-[e4:FRAG]-(f:Mol)
            WHERE e1.prop_synthon=$synthon
            AND EXISTS(e4.prop_pharmfp)
            AND NOT e2.prop_synthon=$first_synthon
            AND NOT e2.prop_synthon=e3.prop_synthon
            WITH usersimilarity.tanimoto_similarity(e1.prop_pharmfp,e4.prop_pharmfp) as sim, e4.prop_synthon as syn, f.smiles as smi, f.cmpd_ids as ids, d.smiles as intermed1, c.smiles as intermed2
            WHERE sim >= $threshold
            AND NOT e4.prop_synthon=$synthon
            RETURN sim, syn, smi, ids, intermed1, intermed2 LIMIT $rev_limit
            """
            % {"smiles": smiles,
               "synthon": synthon,
               "first_synthon": first_synthon,
               "threshold": threshold,
               "rev_limit": rev_limit}
        )

    expansions = []
    intermed_nodes1 = []
    intermed_nodes2 = []
    rev_synthons = []
    rev_similarities = []
    rev_cmpd_ids = []

    with driver.session() as session:
        if strict:
            results = session.run(query, smiles=smiles, synthon=synthon, threshold=threshold, rev_limit=rev_limit)
        else:
            results = session.run(query, smiles=smiles, synthon=synthon, first_synthon=first_synthon, threshold=threshold, rev_limit=rev_limit)
        for res in results:
            expansions.append(res['smi'])
            rev_synthons.append(res['syn'])
            rev_similarities.append(res['sim'])
            rev_cmpd_ids.append(res['ids'])
            intermed_nodes1.append(res['intermed1'])
            if not strict:
                intermed_nodes2.append(res['intermed2'])
    if len(expansions) > 0:
        print('found', len(expansions), 'replacements for', smiles, synthon)
    return expansions, rev_synthons, rev_similarities, rev_cmpd_ids, intermed_nodes1, intermed_nodes2
