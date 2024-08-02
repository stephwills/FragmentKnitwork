import os
from argparse import ArgumentParser

from pymol import cmd
from rdkit import Chem
from tqdm import tqdm


def process_files(dir):
    fragments = os.listdir(dir)
    fragment_dirs = [os.path.join(dir, f) for f in fragments]
    for fragment_dir, frag in tqdm(zip(fragment_dirs, fragments)):
        create_files(fragment_dir, frag)


def create_files(dir, fragment):
    sdf_fname = os.path.join(dir, f"{fragment}.sdf")
    if os.path.exists(sdf_fname):
        mol = Chem.SDMolSupplier(sdf_fname)[0]

    mol_fname = os.path.join(dir, f"{fragment}.mol")
    if not os.path.exists(mol_fname):
        Chem.MolToMolFile(mol, mol_fname)
    else:
        mol = Chem.MolFromMolFile(mol_fname)

    smiles_fname = os.path.join(dir, f"{fragment}_smiles.txt")
    if not os.path.exists(smiles_fname):
        try:
            Chem.RemoveStereochemistry(mol)
            smiles = Chem.MolToSmiles(mol)
            with open(smiles_fname, 'w') as f:
                f.write(smiles)
        except:
            print('No smiles file for', fragment)

    bound_fname = os.path.join(dir, f"{fragment}_bound.pdb")
    apo_fname = os.path.join(dir, f"{fragment}_apo.pdb")
    make_apo_file(bound_fname, apo_fname, False)

    apo_desolv_fname = os.path.join(dir, f"{fragment}_apo-desolv.pdb")
    make_apo_file(bound_fname, apo_desolv_fname, True)


def make_apo_file(pdb_file, new_filename, desolv=False):
    if not os.path.exists(new_filename):
        cmd.reinitialize()
        cmd.load(pdb_file, 'complex')
        cmd.remove('resn dms')
        cmd.remove('resn lig')
        if desolv:
            cmd.remove('solvent')
            cmd.remove('resn hoh')
        cmd.save(new_filename, 'all')
    else:
        print("apo file exists", new_filename)
    return new_filename


if __name__ == "__main__":
    """
    Converts new Fragalysis file format to the old one. Download from Fragalysis selecting the 'save sdf in separate
    subdirectories' option. Run main to get all the necessary files to run the code.
    """
    parser = ArgumentParser()
    parser.add_argument('-d', '--fragalysis_aligned_dir', help="path to aligned dir in fragalysis download")
    args = parser.parse_args()

    process_files(args.fragalysis_aligned_dir)
    print('Files have been processed')
