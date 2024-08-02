
import os
import shutil
import unittest

import FragmentKnitwork.utils.Config as config
from FragmentKnitwork.Fragment.runEnumeration import runEnumeration
from FragmentKnitwork.utils.utils import get_mol, get_protein, get_smiles


class TestFragment(unittest.TestCase):

    def test_run_enumeration(self):
        """
        Test run enumeration runs and writes files to output directory

        :return:
        """
        fA = 'x0310_0A'
        fB = 'x0556_0A'
        fragment_names = [fA, fB]
        target = 'A71EV2A'
        dir = os.path.join(config.WORKING_DIR, 'test/data/Fragalysis/')
        output_dir = os.path.join(config.WORKING_DIR, 'test/data/output')
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        mol_files = [get_mol(target, name, fragalysis_dir=dir) for name in fragment_names]
        pdb_files = [get_protein(target, name, protonated=True, desolv=True, fragalysis_dir=dir) for name in fragment_names]
        fragment_smiles = [get_smiles(target, name, fragalysis_dir=dir) for name in fragment_names]

        runEnumeration(fragment_names, fragment_smiles, target, mol_files, pdb_files, output_dir, ignore_pairs=False,
                       write_smiles_to_file=False, r_group_expansions=True, record_equiv_synthon=True)

        output_files = [i for i in os.listdir(output_dir)]
        tot_files = len(output_files)
        for file in output_files:
            os.remove(os.path.join(output_dir, file))
        shutil.rmtree(output_dir)
        self.assertTrue(tot_files > 0)



if __name__ == '__main__':
    unittest.main()