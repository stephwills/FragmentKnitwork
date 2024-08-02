
import os
import shutil
import unittest

import FragmentKnitwork.utils.Config as config
from FragmentKnitwork.Knitwork.runKnitting import runKnitting
from FragmentKnitwork.Knitwork.reverseQueries import runReverseMerging
from FragmentKnitwork.utils.utils import load_json


class TestKnitwork(unittest.TestCase):

    def test_run_knitting_pure(self):
        """

        :return:
        """
        substructure_dir = os.path.join(config.WORKING_DIR, 'test/data/enumeration')
        substructure_pair_file = os.path.join(substructure_dir, 'substructure_pairs.json')
        target = 'A71EV2A'

        working_dir = os.path.join(config.WORKING_DIR, 'test/data/working')
        output_dir = os.path.join(config.WORKING_DIR, 'test/data/output')

        if not os.path.exists(working_dir):
            os.mkdir(working_dir)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        runKnitting(substructure_pair_file, 1, target, working_dir, output_dir, 10, substructure_dir, 'prop_pharmfp',
                    pure_search=True, prolif_prioritization=False)

        output_file = os.path.join(output_dir, 'x0310_0A-x0556_0A_pure_merge.json')
        data = load_json(output_file)
        len_data = len(data)
        shutil.rmtree(output_dir)
        shutil.rmtree(working_dir)
        self.assertTrue(len_data > 0)


    def test_run_knitting_impure(self):
        """

        :return:
        """
        substructure_dir = os.path.join(config.WORKING_DIR, 'test/data/enumeration')
        substructure_pair_file = os.path.join(substructure_dir, 'substructure_pairs.json')
        target = 'A71EV2A'

        working_dir = os.path.join(config.WORKING_DIR, 'test/data/working')
        output_dir = os.path.join(config.WORKING_DIR, 'test/data/output')

        if not os.path.exists(working_dir):
            os.mkdir(working_dir)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        runKnitting(substructure_pair_file, 1, target, working_dir, output_dir, 10, substructure_dir, 'prop_pharmfp',
                    pure_search=False, prolif_prioritization=False)

        output_file = os.path.join(output_dir, 'x0310_0A-x0556_0A_prop_pharmfp_impure_merge.json')
        data = load_json(output_file)
        len_data = len(data)
        shutil.rmtree(output_dir)
        shutil.rmtree(working_dir)
        self.assertTrue(len_data > 0)


    def test_run_knitting_w_prioritisation(self):
        """

        :return:
        """
        substructure_dir = os.path.join(config.WORKING_DIR, 'test/data/enumeration')
        substructure_pair_file = os.path.join(substructure_dir, 'substructure_pairs.json')
        target = 'A71EV2A'

        working_dir = os.path.join(config.WORKING_DIR, 'test/data/working')
        output_dir = os.path.join(config.WORKING_DIR, 'test/data/output')

        if not os.path.exists(working_dir):
            os.mkdir(working_dir)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        runKnitting(substructure_pair_file, 1, target, working_dir, output_dir, 10, substructure_dir, 'prop_pharmfp',
                    pure_search=False, prolif_prioritization=True, max_prioritize=5)

        output_file = os.path.join(output_dir, 'prolif_prioritized', 'output', 'x0310_0A-x0556_0A_prop_pharmfp_impure_merge.json')
        data = load_json(output_file)
        len_exps = len(data['N#CCC(N)=O*[Xe]c1cc[nH]n1']['expansions'])
        shutil.rmtree(output_dir)
        shutil.rmtree(working_dir)
        self.assertTrue(len_exps == 5)

    def test_run_knitting_w_prioritisation_and_r_group_expansion(self):
        """

        :return:
        """
        substructure_dir = os.path.join(config.WORKING_DIR, 'test/data/enumeration')
        substructure_pair_file = os.path.join(substructure_dir, 'substructure_pairs.json')
        r_group_data_file = os.path.join(substructure_dir, 'r_group_expansions.json')
        equiv_synthon_data_file = os.path.join(substructure_dir, 'equivalent_subnodes.json')
        target = 'A71EV2A'

        working_dir = os.path.join(config.WORKING_DIR, 'test/data/working')
        output_dir = os.path.join(config.WORKING_DIR, 'test/data/output')

        if not os.path.exists(working_dir):
            os.mkdir(working_dir)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        runKnitting(substructure_pair_file, 1, target, working_dir, output_dir, 100, substructure_dir, 'prop_pharmfp',
                    pure_search=False, prolif_prioritization=True, max_prioritize=50, r_group_search=True, r_group_data_file=r_group_data_file,
                    equiv_synthon_data_file=equiv_synthon_data_file)

        output_file = os.path.join(output_dir, 'r_group_expanded', 'x0310_0A-x0556_0A_prop_pharmfp_impure_merge.json')
        data = load_json(output_file)
        len_exps = len(data['N#CCC(N)=O*[Xe]c1cc[nH]n1']['expansions'])
        len_rgroups = len(data['N#CCC(N)=O*[Xe]c1cc[nH]n1']['r_group'])
        shutil.rmtree(output_dir)
        shutil.rmtree(working_dir)
        self.assertTrue(len_exps == len_rgroups)

    def test_reverse_merging_strict(self):
        """

        :return:
        """
        substructure_dir = os.path.join(config.WORKING_DIR, 'test/data/enumeration')
        sdf_file = os.path.join(config.WORKING_DIR, 'test/data/files/impure_round1.sdf')

        working_dir = os.path.join(config.WORKING_DIR, 'test/data/working')
        output_dir = os.path.join(config.WORKING_DIR, 'test/data/output')

        if not os.path.exists(working_dir):
            os.mkdir(working_dir)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        equiv_synthon_data_file = os.path.join(substructure_dir, 'equivalent_subnodes.json')
        runReverseMerging(sdf_file, equiv_synthon_data_file, max_run=2, n_parallel=1, output_dir=output_dir, working_dir=working_dir,
                          threshold=0.9, limit_results=20, strict_linker_search=True)

        output_file = os.path.join(output_dir, 'x0310_0A-x0416_0A_reverse_merge.json')
        data = load_json(output_file)
        shutil.rmtree(output_dir)
        shutil.rmtree(working_dir)
        self.assertTrue(len(data['N#CCC(=O)N[Xe]*[Xe]C1CCCO1*N#CCC(=O)Nc1ncc(C(=O)C2CCCOC2)s1']['expansions']) > 0)



    def test_reverse_merging_loose(self):
        """

        :return:
        """
        substructure_dir = os.path.join(config.WORKING_DIR, 'test/data/enumeration')
        sdf_file = os.path.join(config.WORKING_DIR, 'test/data/files/impure_round1.sdf')

        working_dir = os.path.join(config.WORKING_DIR, 'test/data/working')
        output_dir = os.path.join(config.WORKING_DIR, 'test/data/output')

        if not os.path.exists(working_dir):
            os.mkdir(working_dir)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        equiv_synthon_data_file = os.path.join(substructure_dir, 'equivalent_subnodes.json')
        runReverseMerging(sdf_file, equiv_synthon_data_file, max_run=2, n_parallel=1, output_dir=output_dir,
                          working_dir=working_dir,
                          threshold=0.9, limit_results=20, strict_linker_search=False)

        output_file = os.path.join(output_dir, 'x0310_0A-x0416_0A_reverse_merge.json')
        data = load_json(output_file)
        shutil.rmtree(output_dir)
        shutil.rmtree(working_dir)
        self.assertTrue(len(data['N#CCC(=O)N[Xe]*[Xe]C1CCCO1*N#CCC(=O)Nc1ncc(C(=O)C2CCCOC2)s1']['expansions']) > 0)


if __name__ == '__main__':
    unittest.main()
