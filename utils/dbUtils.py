from FragmentKnitwork.utils.dbConfig import driver


def get_subnodes(smiles: str, terminal_subnodes: bool = False) -> list:
    """
    Get subnodes for a given node (retrieve using SMILES)

    :param smiles: SMILES string for node to retrieve subnodes
    :param terminal_subnodes: whether to only return 'terminal' subnodes (can't be broken down further)
    :return: list of unique subnode SMILES
    """
    terminal_subnode_query = (
        """
        MATCH (a:F2 {smiles: $smiles})-[e:FRAG*0..20]->(f:F2)
        WHERE NOT ()-[:FRAG]-(f)-[:FRAG]->()
        RETURN f
        """
    )
    subnode_query = (
        "MATCH (fa:F2 {smiles: $smiles})-[e:FRAG*0..20]->(f:F2) RETURN f"
    )
    if terminal_subnodes:
        query = terminal_subnode_query
    else:
        query = subnode_query
    with driver.session() as session:
        subnodes = [record['f']['smiles'] for record in session.run(query, smiles=smiles)]
    return list(set(subnodes))


def get_R_groups(smiles: str):
    R_group_query = (
        """
        MATCH (a:F2 {smiles: $smiles})-[e:FRAG*0..20]->(b:F2)
        WHERE NOT ()-[:FRAG]-(b)-[:FRAG]->()
        AND e[-1].prop_synthon contains '[Xe]'
        AND NOT e[-1].prop_synthon=e[-2].prop_synthon
        RETURN e[-1].prop_synthon as synthon, e[-2].prop_synthon as r_group;
        """
    )
    synthons = []
    r_groups = []

    with driver.session() as session:
        results = session.run(R_group_query, smiles=smiles)
        for res in results:
            synthon = res['synthon']
            r_group = res['r_group']
            synthons.append(synthon)
            r_groups.append(r_group)

    return synthons, r_groups


def get_synthons(smiles: str, terminal_synthons: bool = True) -> list:
    """
    Get constituent synthons (compounds added or removed during transformation) for a given node SMILES.
    [Xe] denotes the attachment point.

    :param smiles: SMILES string of node to retrieve synthons
    :param terminal_synthons: whether to return 'terminal' synthons, i.e. can't be broken down more
    :return: list of constituent synthon SMILES strings
    """
    terminal_synthon_query = (
        """
        MATCH (a:F2 {smiles: $smiles})-[e:FRAG*0..20]->(b:F2)
        WHERE NOT ()-[:FRAG]-(b)-[:FRAG]->()
        RETURN e[-1] as edge
        """
    )
    synthon_query = (
        """
        MATCH (a:F2 {smiles: $smiles})-[e:FRAG*0..20]->(b:F2)
        RETURN e[-1] as edge
        """
    )
    if terminal_synthons:
        query = terminal_synthon_query
    else:
        query = synthon_query

    synthons = []
    with driver.session() as session:
        edges = [record["edge"] for record in session.run(query, smiles=smiles)]
        if len(edges) == 0:
            return []
        edges = [edge for edge in edges if edge]
        syns = [edge["prop_synthon"] for edge in edges]
        cores = [edge["prop_core"] for edge in edges]
        synthons.extend(syns)
        synthons.extend(cores)

    synthons = [i for i in synthons if i.count('Xe') == 1]
    return list(set(synthons))
