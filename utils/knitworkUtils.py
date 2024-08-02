from FragmentKnitwork.utils import knitworkConfig as config
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem.Pharm2D.SigFactory import SigFactory


def load_sigFactory(fdef_file=config.FINGERPRINT_FDEF, max_point_count=config.FINGERPRINT_MAXPOINTCOUNT,
                    bins=config.FINGERPRINT_BINS):
    """
    Load signature factory for pharmacophore fp calculation

    :param fdef_file:
    :param max_point_count:
    :param bins:
    :return:
    """
    featFactory = ChemicalFeatures.BuildFeatureFactory(fdef_file)
    sigFactory = SigFactory(featFactory, maxPointCount=max_point_count)
    sigFactory.SetBins(bins)
    sigFactory.Init()
    sigFactory.GetSigSize()
    return sigFactory


def calc_pharm_fp(mol, sigFactory, asStr=True):
    """
    Calculate pharmacophore fingerprint using RDKit

    :param mol:
    :param sigFactory:
    :param asStr:
    :return:
    """
    fp = Generate.Gen2DFingerprint(mol, sigFactory)
    if asStr:
        fp = list(fp)
        return ";".join(map(str, fp))
    else:
        return list(fp)


def get_property_dict(mols):
    """
    Record all the properties of a list of molecules (using names as keys, assuming they are all named).
    This function is used because it seems to be problematic preserving mol properties when running joblib?

    :param mols: list of mols with identical property names
    :return: nested dictionary, {mol_name: {property: val, property: val}, mol_name: {...}, ...}
    """
    property_dict = {}
    names = [mol.GetProp('_Name') for mol in mols]
    props = list(mols[0].GetPropNames())
    # props.remove('_Name')
    for mol, name in zip(mols, names):
        mol_dict = {}
        for prop in props:
            mol_dict[prop] = mol.GetProp(prop)
        property_dict[name] = mol_dict
    return property_dict
