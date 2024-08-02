
import os
import shutil
import unittest

import FragmentKnitwork.utils.Config as config
from FragmentKnitwork.utils.utils import load_json
from FragmentKnitwork.Quilter.quilterPipeline import Pipeline
from FragmentKnitwork.Quilter.quilterWictorPipeline import WictorPipeline


class TestQuilter(unittest.TestCase):

    def test_pipeline(self):
        fA = 'x0310_0A'
        fB = 'x0556_0A'
        target = 'A71EV2A'
        pdb_file = os.path.join(config.WORKING_DIR, 'test/data/Fragalysis/A71EV2A/aligned/A71EV2A-x0310_0A/A71EV2A-x0310_0A_apo-desolv.pdb')
        working_dir = os.path.join(config.WORKING_DIR, 'test/data/working')
        output_dir = os.path.join(config.WORKING_DIR, 'test/data/output')
        scoring_dir = os.path.join(config.WORKING_DIR, 'test/data/scored')
        substructure_dir = os.path.join(config.WORKING_DIR, 'test/data/enumeration')
        file = os.path.join(config.WORKING_DIR, 'test/data/files/x0310_0A-x0556_0A_prop_pharmfp_impure_merge.json')

        if not os.path.exists(working_dir):
            os.mkdir(working_dir)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        if not os.path.exists(scoring_dir):
            os.mkdir(scoring_dir)

            # read in input data

        data = load_json(file)
        smiles, cmpd_ids = [], []
        ref_subnodes, ref_synthons = [], []

        # only used for impure merges & reverse merges
        similarities, used_subnodes, used_synthons = [], [], []

        # read data from input file
        for sub_pair in data:
            print('Substructure pair', sub_pair)
            n_expans = len(data[sub_pair]['expansions'])
            print('N expansions', n_expans)

            # shared between all merge types
            smiles.extend(data[sub_pair]['expansions'])
            cmpd_ids.extend(data[sub_pair]['cmpd_ids'])

            # get ref substructures (not for reverse merge)
            subnode, synthon = sub_pair.split('*')[0], sub_pair.split('*')[1]
            ref_subnodes.extend([subnode] * n_expans)
            ref_synthons.extend([synthon] * n_expans)

            similarities.extend(data[sub_pair]['similarities'])
            used_subnodes.extend([subnode] * n_expans)
            used_synthons.extend(data[sub_pair]['used_synthons'])

        pipeline = Pipeline(smiles,
                            fA,
                            fB,
                            target,
                            pdb_file,
                            cmpd_ids,
                            ref_subnodes,
                            ref_synthons,
                            output_dir,
                            working_dir,
                            scoring_dir,
                            'impure_merge',
                            substructure_dir,
                            similarities=similarities,
                            used_subnodes=used_subnodes,
                            used_synthons=used_synthons,
                            n_cpus=1,
                            move_files=True)
        pipeline.run_filtering()
        pipeline.run_scoring()

        scored_file = os.path.join(scoring_dir, 'scored_data.json')
        data = load_json(scored_file)
        shutil.rmtree(output_dir)
        shutil.rmtree(working_dir)
        shutil.rmtree(scoring_dir)
        self.assertTrue(len(data["names"]) > 0)

    def test_wictor_pipeline(self):
        fA = 'x0310_0A'
        fB = 'x0556_0A'
        target = 'A71EV2A'
        pdb_file = os.path.join(config.WORKING_DIR, 'test/data/Fragalysis/A71EV2A/aligned/A71EV2A-x0310_0A/A71EV2A-x0310_0A_apo-desolv.pdb')
        working_dir = os.path.join(config.WORKING_DIR, 'test/data/working')
        output_dir = os.path.join(config.WORKING_DIR, 'test/data/output')
        scoring_dir = os.path.join(config.WORKING_DIR, 'test/data/scored')
        substructure_dir = os.path.join(config.WORKING_DIR, 'test/data/enumeration')
        file = os.path.join(config.WORKING_DIR, 'test/data/files/x0310_0A-x0556_0A_prop_pharmfp_impure_merge.json')

        if not os.path.exists(working_dir):
            os.mkdir(working_dir)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        if not os.path.exists(scoring_dir):
            os.mkdir(scoring_dir)

            # read in input data

        data = load_json(file)
        smiles, cmpd_ids = [], []
        ref_subnodes, ref_synthons = [], []

        # only used for impure merges & reverse merges
        similarities, used_subnodes, used_synthons = [], [], []

        # read data from input file
        for sub_pair in data:
            print('Substructure pair', sub_pair)
            n_expans = len(data[sub_pair]['expansions'])
            print('N expansions', n_expans)

            # shared between all merge types
            smiles.extend(data[sub_pair]['expansions'])
            cmpd_ids.extend(data[sub_pair]['cmpd_ids'])

            # get ref substructures (not for reverse merge)
            subnode, synthon = sub_pair.split('*')[0], sub_pair.split('*')[1]
            ref_subnodes.extend([subnode] * n_expans)
            ref_synthons.extend([synthon] * n_expans)

            similarities.extend(data[sub_pair]['similarities'])
            used_subnodes.extend([subnode] * n_expans)
            used_synthons.extend(data[sub_pair]['used_synthons'])

        pipeline = WictorPipeline(smiles,
                                  fA,
                                  fB,
                                  target,
                                  pdb_file,
                                  cmpd_ids,
                                  ref_subnodes,
                                  ref_synthons,
                                  output_dir,
                                  working_dir,
                                  scoring_dir,
                                  'impure_merge',
                                  substructure_dir,
                                  similarities=similarities,
                                  used_subnodes=used_subnodes,
                                  used_synthons=used_synthons,
                                  n_cpus=1,
                                  move_files=True,
                                  limit_num_minimize=1)

        pipeline.run_filtering()
        pipeline.run_scoring()

        scored_file = os.path.join(scoring_dir, 'scored_data.json')
        data = load_json(scored_file)
        shutil.rmtree(output_dir)
        shutil.rmtree(working_dir)
        shutil.rmtree(scoring_dir)
        self.assertTrue(len(data["names"]) > 0)


if __name__ == '__main__':
    unittest.main()
