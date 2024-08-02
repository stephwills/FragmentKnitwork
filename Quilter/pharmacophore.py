"""Adapts some code from https://github.com/rdkit/UGM_2016/blob/master/Notebooks/Stiefl_RDKitPh4FullPublication.ipynb"""
import random

from FragmentKnitwork.utils import quilterConfig as config
from FragmentKnitwork.utils.quilterUtils import (get_all_length_combinations,
                                                 get_best,
                                                 get_scoring_function,
                                                 load_feat_factory, loop_list,
                                                 mean_SuCOS)
from FragmentKnitwork.utils.utils import unnest_list
from rdkit import Chem, Geometry
from rdkit.Chem import Mol, rdDistGeom, rdMolTransforms
from rdkit.Chem.Pharm3D import EmbedLib, Pharmacophore
from rdkit.Numerics import rdAlignment

featFactory = load_feat_factory(config.ALIGNMENT_FDEF)


def apply_radii_to_bounds(ph4, radii):
    """

    :param ph4:
    :param radii:
    :return:
    """
    for i in range(len(radii)):
        for j in range(i + 1, len(radii)):
            sumRadii = radii[i] + radii[j]
            ph4.setLowerBound(i, j, max(ph4.getLowerBound(i, j) - sumRadii, 0))
            ph4.setUpperBound(i, j, ph4.getUpperBound(i, j) + sumRadii)


def check_carbons(atom_idxs, mol):
    """
    Check if everything a carbon

    :param atom_idxs:
    :param mol:
    :return:
    """
    atomicNumbers = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in atom_idxs]
    uniqAtomicNumbers = list(set(atomicNumbers))
    if len(uniqAtomicNumbers) == 1 and uniqAtomicNumbers[0] == 6:
        return True
    else:
        return False


def carbon_ring_mappings(feats, ref_substruct, repl_substruct, radii=config.PH4_RADII):
    """
    Generate the possible ring to ring mappings if it is a carbon ring

    :param feats:
    :param ref_substruct:
    :param repl_substruct:
    :param radii:
    :return:
    """
    ref_carbons = check_carbons(feats[0].GetAtomIds(), ref_substruct)  # record if ref ring is all carbons
    pharmacophore = Pharmacophore.Pharmacophore(feats)
    rad = [radii] * len(pharmacophore.getFeatures())
    apply_radii_to_bounds(pharmacophore, rad)

    canMatch, allMatches = EmbedLib.MatchPharmacophoreToMol(repl_substruct,
                                                            featFactory,
                                                            pharmacophore)
    ref_ring_atoms = feats[0].GetAtomIds()
    try:
        if canMatch:
            repl_ring_atoms = []

            for matches in allMatches:
                ringMatch = matches[0]
                ringMatch_carbons = check_carbons(ringMatch.GetAtomIds(),
                                                  repl_substruct)  # check if repl ring is all carbons
                if ref_carbons + ringMatch_carbons >= 1:
                    repl_ring_atoms.append(ringMatch.GetAtomIds())

        all_mappings = [loop_list(atoms) for atoms in repl_ring_atoms]
        all_mappings = unnest_list(all_mappings)

        ring_mappings = [{k: v for k, v in zip(ref_ring_atoms, lst)} for lst in all_mappings]
        return ring_mappings
    except:
        return None


def get_pharmacophore_features_from_ref_mol(ref_mol: Mol, featFactory):
    """
    Get the pharmacophore features from a reference molecule (specifically extracted from the relevant substructure
    used to make the merge).

    :param ref_mol:
    :param featFactory:
    :return:
    """
    pharmacophore_feats = featFactory.GetFeaturesForMol(ref_mol)
    return pharmacophore_feats


def get_transform_matrix(alignRef, confEmbed, atomMatch):
    """

    :param alignRef:
    :param confEmbed:
    :param atomMatch:
    :return:
    """
    alignProbe = []
    for matchIds in atomMatch:
        dummyPoint = Geometry.Point3D(0.0, 0.0, 0.0)
        for id in matchIds:
            dummyPoint += confEmbed.GetAtomPosition(id)
        dummyPoint /= len(matchIds)
        alignProbe.append(dummyPoint)
    return rdAlignment.GetAlignmentTransform(alignRef, alignProbe)


def transform_embeddings(pcophore, embeddings, atom_match, ref_mol, return_SSDs=False, scoring_function=config.ALIGNMENT_SCORE, isFragments=False):
    """
    Borrows code from https://github.com/rdkit/UGM_2016/blob/master/Notebooks/Stiefl_RDKitPh4FullPublication.ipynb.

    :param pcophore:
    :param embeddings:
    :param atom_match:
    :param ref_mol:
    :param return_SSDs:
    :param scoring_function:
    :param isFragments:
    :return:
    """
    func = get_scoring_function(scoring_function)
    alignRef = [f.GetPos() for f in pcophore.getFeatures()]
    SSDs = []
    scores = []
    for embedding in embeddings:
        conf = embedding.GetConformer()
        SSD, transformMatrix = get_transform_matrix(alignRef, conf, atom_match)
        rdMolTransforms.TransformConformer(conf, transformMatrix)
        SSDs.append(SSD)
        if not isFragments:
            scores.append(func(ref_mol, embedding))
        if isFragments:  # if ref_mol is (fragmentA, fragmentB) and we want to check sucos against them
            score = mean_SuCOS(embedding, ref_mol[0], ref_mol[1], asMol=True, returnFragmentMols=False)
            scores.append(score)

    if return_SSDs:
        return SSDs, scores
    else:
        return scores


def generate_embeddings(pharmacophore, to_embed, radii_dist, embeddings_count=config.EMBEDDINGS_COUNT):
    """
    Attempt to generate embedding of a molecule using a Pharmacophore.Pharmacophore object. Adapts code from
    https://github.com/rdkit/UGM_2016/blob/master/Notebooks/Stiefl_RDKitPh4FullPublication.ipynb.

    :param pharmacophore: Pharmacophore object created from feats
    :param to_embed: the molecule to embed
    :param radii_dist: radii_dist
    :param embeddings_count: number of embeddings to try
    :return:
    """
    radii = [radii_dist] * len(pharmacophore.getFeatures())
    apply_radii_to_bounds(pharmacophore, radii)
    canMatch, allMatches = EmbedLib.MatchPharmacophoreToMol(Chem.Mol(to_embed),
                                                            featFactory,
                                                            pharmacophore)
    if not canMatch:
        return None, None, None

    boundsMat = rdDistGeom.GetMoleculeBoundsMatrix(to_embed)
    (failed, boundsMatMatches, matched, matchDetails) = EmbedLib.MatchPharmacophore(allMatches,
                                                                                    boundsMat,
                                                                                    pharmacophore,
                                                                                    useDownsampling=True)

    if failed == 0:  # if embedding was successful
        atomMatch = [list(x.GetAtomIds()) for x in matched]
        to_embed = Chem.AddHs(to_embed, addCoords=True)
        bm, embeddings, numFail = EmbedLib.EmbedPharmacophore(to_embed,
                                                              atomMatch,
                                                              pharmacophore,
                                                              count=embeddings_count)

        if len(embeddings) > 0:
            return embeddings, matched, atomMatch

    return None, None, None


def embed_substructure_with_pharmacophores(ref_substruct: Mol, used_substruct: Mol, radii_dist=config.PH4_RADII,
                                           scoring_mode=config.SCORING_MODE):
    """
    Embed a replacement substructure using the pharmacophore features from the original substructure. The exact
    embedding and combination of pharmacophore features chosen according to best-scoring (right now using SuCOS).

    :param ref_substruct:
    :param used_substruct:
    :param radii_dist:
    :param scoring_mode:
    :return:
    """
    # get all possible combinations of pharmacophore features from the reference substructure
    pharmacophore_feats = get_pharmacophore_features_from_ref_mol(ref_substruct, featFactory)
    combinations_feats = get_all_length_combinations(pharmacophore_feats)

    all_embeddings = []
    all_values = []
    all_ref_atom_ids = []
    all_repl_atom_ids = []

    all_feats = []

    # check for all-carbon rings
    rings = ['Aliphatic', 'Aromatic']
    for feats in combinations_feats:
        pharmacophore = Pharmacophore.Pharmacophore(feats)
        embeddings, matched, atomMatch = generate_embeddings(pharmacophore, used_substruct, radii_dist)

        if embeddings:
            # now that we have a successful embedding - record the atom matches
            ref_atom_ids = [f.GetAtomIds() for f in feats]
            repl_atom_ids = [m.GetAtomIds() for m in matched]
            all_ref_atom_ids.append(ref_atom_ids)
            all_repl_atom_ids.append(repl_atom_ids)
            scores = transform_embeddings(pharmacophore, embeddings, atomMatch, ref_substruct)
            best_score, idx, best_embedding = get_best(embeddings, scores, scoring_mode)
            all_embeddings.append(best_embedding)
            all_values.append(best_score)
            all_feats.append(feats)

    if len(all_embeddings) == 0:
        return None, None, None, None

    best_value, best_idx, best_embedding = get_best(all_embeddings, all_values, scoring_mode)
    best_feats = all_feats[best_idx]

    # get the atom ids of the features used to get the best embedding
    best_ref_atom_ids = all_ref_atom_ids[best_idx]
    best_repl_atom_ids = all_repl_atom_ids[best_idx]

    # a bit complicated but if a carbon ring, there should be multiple possible atom-to-atom mappings
    # this generates the possible mappings
    if len(best_feats) == 1 and best_feats[0].GetFamily() in rings:  # check if only one feature used and if it is a ring ph4
        carbon_ring_maps = carbon_ring_mappings(best_feats, ref_substruct, used_substruct)
        if carbon_ring_maps:
            if len(carbon_ring_maps) > 0:
                return best_embedding, best_ref_atom_ids, best_repl_atom_ids, carbon_ring_maps
            else:
                return best_embedding, best_ref_atom_ids, best_repl_atom_ids, False
        else:
            return best_embedding, best_ref_atom_ids, best_repl_atom_ids, False

    else:
        return best_embedding, best_ref_atom_ids, best_repl_atom_ids, False


def embed_mol_with_pharmacophores(fragmentA: Mol, fragmentB: Mol, merge: Mol, radii_dist: float = config.PH4_RADII,
                                  scoring_mode=config.SCORING_MODE, n_feat_combinations: int = config.N_FEAT_COMBINATIONS):
    """
    Embed a molecule using the pharmacophores from the two inspiration fragments

    :param fragmentA:
    :param fragmentB:
    :param merge:
    :param radii_dist:
    :param scoring_mode:
    :param n_feat_combinations:
    :return:
    """
    pharmacophore_featsA = get_pharmacophore_features_from_ref_mol(fragmentA, featFactory)
    pharmacophore_featsB = get_pharmacophore_features_from_ref_mol(fragmentB, featFactory)
    # combine pharmacophore features from both fragments
    pharmacophore_feats = pharmacophore_featsA + pharmacophore_featsB
    # get all possible pharmacophore combinations
    combinations_feats = get_all_length_combinations(pharmacophore_feats)

    all_embeddings = []
    all_values = []
    all_ref_atom_ids = []
    all_repl_atom_ids = []
    all_feats = []

    # as it may take too long to sample all possible pharmacophore combinations, take a random selection (default 100)
    random.seed(42)
    combinations_feats = random.choices(combinations_feats, k=n_feat_combinations)

    for feats in combinations_feats:
        # for each combination of pharmacophores
        pharmacophore = Pharmacophore.Pharmacophore(feats)  # create pharmacophore object
        embeddings, matched, atomMatch = generate_embeddings(pharmacophore, merge, radii_dist)

        if embeddings:
            # get the best embedding for each ph4 combinations
            # now that we have a successful embedding - record the atom matches
            ref_atom_ids = [f.GetAtomIds() for f in feats]
            repl_atom_ids = [m.GetAtomIds() for m in matched]
            all_ref_atom_ids.append(ref_atom_ids)
            all_repl_atom_ids.append(repl_atom_ids)
            # get SuCOS scores associated with embeddings
            scores = transform_embeddings(pharmacophore, embeddings, atomMatch, [fragmentA, fragmentB], isFragments=True)
            best_score, idx, best_embedding = get_best(embeddings, scores, scoring_mode)
            all_embeddings.append(best_embedding)
            all_values.append(best_score)
            all_feats.append(feats)

    if len(all_embeddings) == 0:
        return None, None, None, None

    # get the best overall ph4 combination
    best_value, best_idx, best_embedding = get_best(all_embeddings, all_values, scoring_mode)
    best_feats = all_feats[best_idx]

    # get the atom ids of the features used to get the best embedding
    best_ref_atom_ids = all_ref_atom_ids[best_idx]
    best_repl_atom_ids = all_repl_atom_ids[best_idx]

    return best_embedding, best_ref_atom_ids, best_repl_atom_ids, best_feats
