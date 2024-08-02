
import itertools
import logging
import os
import shutil
import time

from fragmenstein.faux_victors import Wictor
from FragmentKnitwork.Quilter.map import get_custom_maps
from FragmentKnitwork.utils import quilterConfig as config
from FragmentKnitwork.utils.quilterUtils import (
    eval_conformer, get_best, run_attach_atom_check,
    run_attach_atom_check_for_loose_expansion)
from FragmentKnitwork.utils.utils import disable_rdlogger, load_json
from rdkit import Chem

disable_rdlogger()


def alignment(smiles, name, merge, ref_subnodes, ref_synthons, used_subnode, used_synthon, pdb_file, output_dir,
              minimize=config.PYROSETTA_MINIMIZATION, scoring_function=config.ALIGNMENT_SCORE, mode=config.SCORING_MODE,
              rmsd_threshold=config.RMSD_THRESHOLD, ddg_filter=config.DDG_FILTER, remove_wictor=config.REMOVE_WICTOR,
              reverse_merge=False, covalent_resi=config.COVALENT_RESI):
    """
    This is the function to perform alignment of bioisosteric merges (whereby one or two substructures have been
    replaced). Starts by using pharmacophores to place just the substructure and using this to provide Fragmenstein
    with a custom mapping. Can perform minimization with PyRosetta through Fragmenstein. Code runs Wictor (no minimization)
    for all possible substructure combinations and then either keeps or runs minimization on the best one.

    :param smiles: SMILES of merge
    :param name: name of merge
    :param merge: merge Mol
    :param ref_subnodes: list of possible ref subnode mols (could be multiple matches)
    :param ref_synthons: list of possible ref synthon mols (could be multiple matches)
    :param used_subnode: SMILES of the actual subnode incorporated (only different if reverse merge)
    :param used_synthon: SMILES of actual synthon incorporated
    :param pdb_file: pdb file of the protein for placement
    :param output_dir: where to save output files
    :param minimize: whether to perform minimization
    :param scoring_function: scoring function to use for scoring (only SuCOS used atm)
    :param mode: scoring mode
    :param rmsd_threshold:
    :param ddg_filter:
    :param remove_wictor:
    :param reverse_merge: whether the merge is a 'reverse merge' (two substructures have been replaced rather than one)
    :return:
    """
    # record time taken
    start = time.time()
    time_taken = None
    print('Running alignment for', name)

    name_hyph = name.replace('_', '-')  # fragmenstein saves with hyph instead of underscore

    try:
        # get all possible combinations of substructure for placement
        used_subnode, used_synthon = Chem.MolFromSmiles(used_subnode), Chem.MolFromSmiles(used_synthon)
        ref_combinations = list(itertools.product(ref_subnodes, ref_synthons))
        wictor_dirs, ref_combs_run, mol_files, scores, custom_maps = [], [], [], [], []

        for p, ref_comb in enumerate(ref_combinations):
            # print(f'Running substructure combination {p + 1} of {len(ref_combinations)}')
            ref_subnode, ref_synthon = ref_comb[0], ref_comb[1]
            custom_maps_for_comb = get_custom_maps(merge, ref_subnode, ref_synthon,
                                                   used_subnode, used_synthon, reverse_merge=reverse_merge)

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

        # take time if no minimization performed
        end = time.time()
        time_taken = round(end-start, 2)

        if len(mol_files) == 0:
            print(name, 'fail')
            return False, None, None, None, time_taken

        _, best_idx, best_map = get_best(custom_maps, scores, mode)
        mapping = custom_maps[best_idx]
        keep_wictor_dir = wictor_dirs[best_idx]
        ref_subnode, ref_synthon = ref_combs_run[best_idx][0], ref_combs_run[best_idx][1]

        if not minimize:
            if remove_wictor:
                wictor_dirs.remove(keep_wictor_dir)
                [shutil.rmtree(dir) for dir in wictor_dirs]
            return True, None, keep_wictor_dir, None, time_taken

        else:
            if remove_wictor:
                [shutil.rmtree(dir) for dir in wictor_dirs]

            # run minimization using Victor
            from fragmenstein import Igor, Victor
            Igor.init_pyrosetta()
            v = Victor(hits=[ref_subnode, ref_synthon], pdb_filename=pdb_file, covalent_resi=covalent_resi)
            # print(f'Running fragmenstein for best map')
            custom_map = mapping
            v.work_path = output_dir
            v.enable_stdout(level=logging.FATAL)
            v.place(smiles, long_name=f"{name}", custom_map=custom_map)
            merge_dir = os.path.join(output_dir, name_hyph)

            # record time if minimization performed
            end = time.time()
            time_taken = round(end - start, 2)

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
        print('Error', e, 'for merge', name)
        return None, None, None, str(e), time_taken


def placement(smiles, name, ref_subnodes, ref_synthons, pdb_file, output_dir, minimize=config.PYROSETTA_MINIMIZATION,
              scoring_function=config.ALIGNMENT_SCORE, mode=config.SCORING_MODE, rmsd_threshold=config.RMSD_THRESHOLD,
              ddg_filter=config.DDG_FILTER, remove_wictor=config.REMOVE_WICTOR, covalent_resi=config.COVALENT_RESI):
    """
    This is the function to perform alignment of pure merges (exact substructures have been incorporated). Places the
    merge SMILES directly using the different combinations of substructure and saves the best-scoring. Can perform
    minimization with PyRosetta through Fragmenstein. Code runs Wictor (no minimization) for all possible substructure
    combinations and then either keeps or runs minimization on the best one.

    :param smiles: SMILES of merge
    :param name: name of merge
    :param ref_subnodes: list of possible ref subnode mols (could be multiple matches)
    :param ref_synthons: list of possible ref synthon mols (could be multiple matches)
    :param pdb_file: pdb file of the protein for placement
    :param output_dir: where to save output files
    :param minimize: whether to perform minimization
    :param scoring_function: scoring function to use for scoring (only SuCOS used atm)
    :param mode: scoring mode
    :param rmsd_threshold:
    :param ddg_filter:
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

            # take time if no minimization performed
            end = time.time()
            time_taken = round(end - start, 2)

            if len(mol_files) == 0:
                print(name, 'fail')
                return False, None, None, None, time_taken

        _, best_idx, keep_wictor_dir = get_best(wictor_dirs, scores, mode)
        ref_subnode, ref_synthon = ref_combs_run[best_idx][0], ref_combs_run[best_idx][1]

        if not minimize:
            if remove_wictor:
                wictor_dirs.remove(keep_wictor_dir)
                [shutil.rmtree(dir) for dir in wictor_dirs]
            return True, None, keep_wictor_dir, None, time_taken

        else:
            if remove_wictor:
                [shutil.rmtree(dir) for dir in wictor_dirs]

            # run minimization using Victor
            from fragmenstein import Igor, Victor
            Igor.init_pyrosetta()
            v = Victor(hits=[ref_subnode, ref_synthon], pdb_filename=pdb_file, covalent_resi=covalent_resi)
            v.work_path = output_dir
            v.enable_stdout(level=logging.FATAL)
            v.place(smiles, long_name=f"{name}")
            # new_name = name.replace('_', '-')
            merge_dir = os.path.join(output_dir, name_hyph)

            # record time if minimization performed
            end = time.time()
            time_taken = round(end - start, 2)

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
        print('Error', e, 'for merge', name)
        return None, None, None, str(e), time_taken


def reverse_alignment(smiles, original_smiles, intermed_smiles1, name, merge, ref_subnodes, ref_synthons,
                      used_subnode, used_synthon, pdb_file, output_dir, intermed_smiles2=None,
                      minimize=config.PYROSETTA_MINIMIZATION, scoring_function=config.ALIGNMENT_SCORE,
                      mode=config.SCORING_MODE, rmsd_threshold=config.RMSD_THRESHOLD, ddg_filter=config.DDG_FILTER,
                      remove_wictor=config.REMOVE_WICTOR, covalent_resi=config.COVALENT_RESI):
    """
    Function for running alignment of a merge whereby both susbstructures have been replaced. First checks that
    the attachment atom is in the right place (i.e. the original substructure in the first bioisosteric merge
    was removed and a replacement added in the same place).

    :param smiles:
    :param original_smiles:
    :param intermed_smiles1:
    :param name:
    :param merge:
    :param ref_subnodes:
    :param ref_synthons:
    :param used_subnode:
    :param used_synthon:
    :param pdb_file:
    :param output_dir:
    :param intermed_smiles2:
    :param minimize:
    :param scoring_function:
    :param mode:
    :param rmsd_threshold:
    :param ddg_filter:
    :param remove_wictor:
    :return:
    """
    original_mol = Chem.MolFromSmiles(original_smiles)
    intermed_mol1 = Chem.MolFromSmiles(intermed_smiles1)
    exp_mol = Chem.MolFromSmiles(smiles)
    if intermed_smiles2:  # this is for non-strict queries (whereby a change in linker occurred)
        intermed_mol2 = Chem.MolFromSmiles(intermed_smiles2)
        res = run_attach_atom_check_for_loose_expansion(intermed_mol2, intermed_mol1, exp_mol)
    else:  # this is for strict queries (whereby there was no change in the linker)
        res = run_attach_atom_check(original_mol, intermed_mol1, exp_mol)
    if not res:
        print('Attachment atom fail')
        return False, None, None, 'attachment atom fail', None
    else:
        return alignment(smiles, name, merge, ref_subnodes, ref_synthons, used_subnode, used_synthon, pdb_file,
                         output_dir, minimize=minimize, scoring_function=scoring_function, mode=mode,
                         rmsd_threshold=rmsd_threshold, ddg_filter=ddg_filter, remove_wictor=remove_wictor,
                         reverse_merge=True, covalent_resi=covalent_resi)


def r_group_alignment(smiles, original_mol, output_dir, pdb_file, name, covalent_resi=config.COVALENT_RESI):
    """
    Code to run placement for an R-group expansion (substituent added to an existing merge), using the original
    placed molecule as constraint for the placement.

    :param smiles: SMILES of the merge after R-group expansion
    :param original_mol: original placed Mol before R-group expansion
    :param output_dir: output directory
    :param pdb_file: pdb file of protein for placement
    :param name: name of the new merge
    :return:
    """
    # record time taken
    start = time.time()

    # use the original placed mol for placement
    original_mol.SetProp('_Name', 'hit')
    name_hyph = name.replace('_', '-')
    merge_dir = os.path.join(output_dir, name_hyph)
    if os.path.exists(merge_dir):
        print('Placement for R-group expansion already run', name)
        return True, Chem.MolFromMolFile(os.path.join(merge_dir, f"{name_hyph}.minimised.mol"))

    # perform placement with minimization using Victor
    from fragmenstein import Igor, Victor
    Igor.init_pyrosetta()
    v = Victor(hits=[original_mol], pdb_filename=pdb_file, covalent_resi=covalent_resi)
    v.work_path = output_dir
    v.enable_stdout(level=logging.FATAL)
    v.place(smiles, long_name=name)

    # record time taken
    end = time.time()
    time_taken = round(end-start, 2)

    # check if successful (and mol file exists)
    mol_file = os.path.join(merge_dir, f"{name_hyph}.minimised.mol")
    if os.path.exists(mol_file):
        print('pass', name_hyph)
        return True, Chem.MolFromMolFile(mol_file), time_taken
    else:
        print('fail', name_hyph)
        return False, None, time_taken
