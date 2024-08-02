
import itertools
import logging
import os
import shutil
import time

from fragmenstein.faux_victors import Wictor
from FragmentKnitwork.utils import quilterConfig as config
from FragmentKnitwork.utils.quilterUtils import eval_conformer, get_best
from FragmentKnitwork.utils.utils import disable_rdlogger, load_json
from FragmentKnitwork.Quilter.map import get_custom_maps
from rdkit import Chem

disable_rdlogger()

# TODO: currently has a lot of duplicated code from align.py

def wictor_alignment(smiles, name, merge, ref_subnodes, ref_synthons, used_subnode, used_synthon, pdb_file, output_dir,
                     scoring_function=config.ALIGNMENT_SCORE, mode=config.SCORING_MODE, remove_wictor=config.REMOVE_WICTOR):
    """
    This is the function to perform initial Victor alignment of bioisosteric merges (whereby one substructure has been
    replaced). Will use this to select the best merges later for minimization.

    :param smiles: SMILES of merge
    :param name: name of merge
    :param merge: merge Mol
    :param ref_subnodes: list of possible ref subnode mols (could be multiple matches)
    :param ref_synthons: list of possible ref synthon mols (could be multiple matches)
    :param used_subnode: SMILES of the actual subnode incorporated (only different if reverse merge)
    :param used_synthon: SMILES of actual synthon incorporated
    :param pdb_file: pdb file of the protein for placement
    :param output_dir: where to save output files
    :param scoring_function: scoring function to use for scoring (only SuCOS used atm)
    :param mode:
    :param remove_wictor:
    :return:
    """
    # record time taken
    start = time.time()
    time_taken = None
    print('Running wictor alignment for', name)

    name_hyph = name.replace('_', '-')

    try:
        # get all possible combinations of substructure for placement
        used_subnode = Chem.MolFromSmiles(used_subnode)
        used_synthon = Chem.MolFromSmiles(used_synthon)
        ref_combinations = list(itertools.product(ref_subnodes, ref_synthons))

        wictor_dirs, ref_combs_run, mol_files, scores, custom_maps = [], [], [], [], []

        for p, ref_comb in enumerate(ref_combinations):
            # print(f'Running substructure combination {p + 1} of {len(ref_combinations)}')
            ref_subnode, ref_synthon = ref_comb[0], ref_comb[1]
            custom_maps_for_comb = get_custom_maps(merge, ref_subnode, ref_synthon, used_subnode, used_synthon,
                                                   reverse_merge=False)

            if custom_maps_for_comb:
                ref_subnode.SetProp('_Name', 'ref_subnode')
                ref_synthon.SetProp('_Name', 'ref_synthon')

                for i, mapping in enumerate(custom_maps_for_comb):
                    map_name = f"{name_hyph}-map-{p}-{i}"
                    wictor_dir = os.path.join(output_dir, map_name)
                    w = Wictor(hits=[ref_subnode, ref_synthon], pdb_filename=pdb_file)
                    # print(f'Running fragmenstein for custom map {i + 1} of {len(custom_maps_for_comb)}')
                    custom_map = {"ref_subnode": mapping[0],
                                  "ref_synthon": mapping[1]}
                    w.work_path = output_dir
                    w.enable_stdout(level=logging.FATAL)
                    try:
                        w.place(smiles, long_name=map_name, custom_map=custom_map)
                        mol_file = os.path.join(output_dir, map_name, f"{map_name}.minimised.mol")
                        scores.append(eval_conformer(Chem.MolFromMolFile(mol_file),
                                                     ref_subnode, ref_synthon, scoring_function))
                        mol_files.append(mol_file)
                        custom_maps.append(custom_map)
                        wictor_dirs.append(wictor_dir)
                        ref_combs_run.append(ref_comb)

                    except Exception as e:
                        if os.path.exists(wictor_dir):
                            shutil.rmtree(wictor_dir)
                        # print('Failed for', map_name)
                        print(e)
            else:
                print('no custom map')

        # record time taken
        end = time.time()
        time_taken = round(end-start, 2)

        if len(mol_files) == 0:
            print(name, 'fail')
            return False, None, None, None, None, time_taken

        # get the best scoring after Wictor placement
        best_score, best_idx, best_map = get_best(custom_maps, scores, mode)  # return best_score, idx, items[idx]
        mapping = custom_maps[best_idx]
        ref_subnode, ref_synthon = ref_combs_run[best_idx][0], ref_combs_run[best_idx][1]

        if remove_wictor:
            [shutil.rmtree(dir) for dir in wictor_dirs]

        return True, best_score, ref_subnode, ref_synthon, mapping, time_taken

    except Exception as e:
        print('Error', e, 'for merge', name)
        return False, None, None, None, None, time_taken


def wictor_placement(smiles, name, ref_subnodes, ref_synthons, pdb_file, output_dir, scoring_function=config.ALIGNMENT_SCORE,
                     mode=config.SCORING_MODE, remove_wictor=config.REMOVE_WICTOR):
    """
    This is the function to perform initial Victor alignment of pure merges (whereby no substructure has been
    replaced). Will use this to select the best merges later for minimization.

    :param smiles: SMILES of merge
    :param name: name of merge
    :param ref_subnodes: list of possible ref subnode mols (could be multiple matches)
    :param ref_synthons: list of possible ref synthon mols (could be multiple matches)
    :param pdb_file: pdb file of the protein for placement
    :param output_dir: where to save output files
    :param scoring_function: scoring function to use for scoring (only SuCOS used atm)
    :param mode:
    :param remove_wictor:
    :return:
    """
    start = time.time()
    time_taken = None
    print('Running placement for', name)

    name_hyph = name.replace('_', '-')

    try:
        # get all possible combinations of substructure for placement
        ref_combinations = list(itertools.product(ref_subnodes, ref_synthons))
        wictor_dirs, ref_combs_run, mol_files, scores = [], [], [], []

        for p, ref_comb in enumerate(ref_combinations):
            # print(f'Running substructure combination {p + 1} of {len(ref_combinations)}')
            ref_subnode, ref_synthon = ref_comb[0], ref_comb[1]
            ref_subnode.SetProp('_Name', 'ref_subnode')
            ref_synthon.SetProp('_Name', 'ref_synthon')

            comb_name = f"{name_hyph}-comb-{p}"
            wictor_dir = os.path.join(output_dir, comb_name)
            w = Wictor(hits=[ref_subnode, ref_synthon], pdb_filename=pdb_file)
            w.work_path = output_dir
            w.enable_stdout(level=logging.FATAL)

            try:
                w.place(smiles, long_name=comb_name)
                mol_file = os.path.join(output_dir, comb_name, f"{comb_name}.minimised.mol")
                scores.append(eval_conformer(Chem.MolFromMolFile(mol_file), ref_subnode, ref_synthon, scoring_function))
                mol_files.append(mol_file)
                wictor_dirs.append(wictor_dir)
                ref_combs_run.append(ref_comb)

            except Exception as e:
                if os.path.exists(wictor_dir):
                    shutil.rmtree(wictor_dir)

            # record time taken
            end = time.time()
            time_taken = round(end-start, 2)

            if len(mol_files) == 0:
                print(name, 'fail')
                return False, None, None, None, None, time_taken

        best_score, best_idx, keep_wictor_dir = get_best(wictor_dirs, scores, mode)
        ref_subnode, ref_synthon = ref_combs_run[best_idx][0], ref_combs_run[best_idx][1]

        if remove_wictor:
            [shutil.rmtree(dir) for dir in wictor_dirs]

        return True, best_score, ref_subnode, ref_synthon, None, time_taken

    except Exception as e:
        print('Error', e, 'for merge', name)
        return None, None, None, None, None, time_taken


def victor_alignment(smiles, name, ref_subnode, ref_synthon, mapping, pdb_file, output_dir, rmsd_threshold=config.RMSD_THRESHOLD,
                     ddg_filter=config.DDG_FILTER, covalent_resi=config.COVALENT_RESI):
    """
    This is the code to perform alignment with minimization after selecting the best Wictor alignments.

    :param smiles:
    :param name:
    :param ref_subnode:
    :param ref_synthon:
    :param mapping:
    :param pdb_file:
    :param output_dir:
    :param rmsd_threshold:
    :param ddg_filter:
    :param: covalent_resi:
    :return:
    """
    start = time.time()
    time_taken = None

    name_hyph = name.replace('_', '-')

    try:
        from fragmenstein import Igor, Victor
        Igor.init_pyrosetta()
        ref_subnode.SetProp('_Name', 'ref_subnode')
        ref_synthon.SetProp('_Name', 'ref_synthon')
        v = Victor(hits=[ref_subnode, ref_synthon], pdb_filename=pdb_file)
        v.covalent_resi=covalent_resi
        v.work_path = output_dir
        v.enable_stdout(level=logging.FATAL)

        if mapping:
            new_custom_map = {"ref_subnode": {}, "ref_synthon": {}}
            for key in mapping["ref_subnode"].keys():
                new_custom_map["ref_subnode"][int(key)] = mapping["ref_subnode"][key]
            for key in mapping["ref_synthon"].keys():
                new_custom_map["ref_synthon"][int(key)] = mapping["ref_synthon"][key]
            v.place(smiles, long_name=f"{name}", custom_map=new_custom_map)

        else:
            v.place(smiles, long_name=f"{name}")

        # record time taken
        end = time.time()
        time_taken = round(end-start, 2)

        merge_dir = os.path.join(output_dir, name_hyph)

        if ddg_filter or rmsd_threshold:
            json_file = os.path.join(merge_dir, f"{name_hyph}.minimised.json")
            data = load_json(json_file)
            deltaG = v.ddG
            comRMSD = data["mRMSD"]  # RMSD between two fragments and merge

            results = []
            if ddg_filter:
                res = deltaG < 0
                results.append(res)
            if rmsd_threshold:
                res = comRMSD <= rmsd_threshold
                results.append(res)

            if len(set(results)) == 1 and results[0]:
                print(name, 'pass')
                return True, v.minimized_mol, merge_dir, None, time_taken
            else:
                shutil.rmtree(merge_dir)
                print(name, 'fail')
                return False, None, None, None, time_taken

        else:
            print(name, 'pass')
            return True, v.minimized_mol, merge_dir, None, time_taken

    except Exception as e:
        print(e)
        return False, None, None, str(e), time_taken
