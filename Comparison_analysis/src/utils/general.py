import pathlib

import pandas
import pandas as pd
from src import definitions
import json
import requests

from src.utils.db import query_activitydb, query_activities
from src.utils.getter import get_activitydb_pdb_ligand_table, get_activitydb_file_processing, get_full_chembl_mol_df, \
    get_chembl_id_to_hetcode_dict


def split_combinations_by_target_found(current_combinations, recalc):
    activityfinder_target_path = (
            definitions.output_files_path / "activityfinder_target.json"
    )
    combinations_to_look_into_target = []
    investigate_combinations_path_target = (
            definitions.output_files_path / "combinations_to_look_into_target2.json"
    )
    combinations_to_look_into_notarget = []
    investigate_combinations_path_notarget = (
            definitions.output_files_path / "combinations_to_look_into_notarget2.json"
    )
    if (
        investigate_combinations_path_target.exists()
        and investigate_combinations_path_notarget.exists()
    ) and not recalc:
        with open(investigate_combinations_path_target, "r") as f:
            combinations_to_look_into_target = json.load(f)
        with open(investigate_combinations_path_notarget, "r") as f:
            combinations_to_look_into_notarget = json.load(f)
    else:
        if not activityfinder_target_path.exists():
            raise RuntimeError("Activityfinder target data has not been generated yet")
        with open(
                definitions.output_files_path / "activityfinder_target.json", "r"
        ) as f:
            activityfinder_targets = json.load(f)
        for missing_combination in current_combinations:
            found = False
            pdb, target, molecule = missing_combination.split(":")
            if pdb in activityfinder_targets:
                if target in activityfinder_targets[pdb]:
                    found = True
            if found:
                combinations_to_look_into_target.append([pdb, target, molecule])
            else:
                combinations_to_look_into_notarget.append([pdb, target, molecule])
        with open(investigate_combinations_path_target, "w") as f:
            json.dump(combinations_to_look_into_target, f)
        with open(investigate_combinations_path_notarget, "w") as f:
            json.dump(combinations_to_look_into_notarget, f)
    print(
        f"Combinations to look into with target, we have {len(combinations_to_look_into_target)} ones to investigate"
    )
    print(
        f"Combinations to look into without target, we have {len(combinations_to_look_into_notarget)} ones to investigate"
    )
    return combinations_to_look_into_target, combinations_to_look_into_notarget


def extract_chembl_inchi_smiles(recalc, combinations_to_investigate):
    chembl_smiles_dict_path = definitions.output_files_path / "chembl_smiles_dict.json"
    chembl_inchi_dict_path = definitions.output_files_path / "chembl_inchi_dict.json"
    if (
        chembl_smiles_dict_path.exists()
        and chembl_inchi_dict_path.exists()
        and not recalc
    ):
        with open(chembl_smiles_dict_path, "r") as f:
            chembl_smiles_dict = json.load(f)
        with open(chembl_inchi_dict_path, "r") as f:
            chembl_inchi_dict = json.load(f)
    else:
        # Calc the full chembl molecule dataframe
        molecule_df = get_full_chembl_mol_df(recalc)
        all_molecules = [x[2] for x in combinations_to_investigate]
        chembl_smiles_dict = {}
        chembl_inchi_dict = {}
        for index, row in molecule_df.iterrows():
            chembl_id = row["chembl_id"]
            standard_inchi = row["standard_inchi"]
            can_smiles = row["canonical_smiles"]
            if chembl_id in all_molecules:
                chembl_smiles_dict[chembl_id] = can_smiles
                chembl_inchi_dict[chembl_id] = standard_inchi
        with open(chembl_smiles_dict_path, "w") as f:
            json.dump(chembl_smiles_dict, f)
        with open(chembl_inchi_dict_path, "w") as f:
            json.dump(chembl_inchi_dict, f)
    print("Chembl smiles and inchi found")
    return chembl_smiles_dict, chembl_inchi_dict


def generate_smiles_for_combinations(combinations_to_investigate, recalc_all=False):
    chembl_smiles_dict, _ = extract_chembl_inchi_smiles(recalc_all, combinations_to_investigate)
    chembl_id_to_hetcode = get_chembl_id_to_hetcode_dict(recalc_all)
    activitydb_liganddata = get_activitydb_pdb_ligand_table(recalc_all)
    hetcode_to_our_smiles = {
        row["ligand_id"]: row["naomi_unique_smiles_chiral"]
        for index, row in activitydb_liganddata.iterrows()
        if row["ligand_id"] in [chembl_id_to_hetcode[x] for x in chembl_id_to_hetcode]
    }

    final_molecule_data = []
    smiles_data = []
    all_hetcodes = []
    hetcode_not_found = 0
    for combination in combinations_to_investigate:
        pdb, target, molecule = combination
        hetcode = chembl_id_to_hetcode[molecule]
        if not hetcode in hetcode_to_our_smiles:
            hetcode_not_found += 1
            our_smiles = ""
            chembl_smiles = ""
        else:
            our_smiles = hetcode_to_our_smiles[hetcode]
            chembl_smiles = chembl_smiles_dict[molecule]
        to_append = [pdb, target, molecule, hetcode, chembl_smiles, our_smiles]
        name = pdb + "-" + target + "-" + molecule + "-" + hetcode
        smiles_data.append([chembl_smiles, name])
        smiles_data.append([our_smiles, name])
        all_hetcodes.append(hetcode)
        final_molecule_data.append(to_append)
    final_molecule_data = pd.DataFrame(
        final_molecule_data,
        columns=[
            "PDB",
            "Chembl Target",
            "Chembl Molecule",
            "hetcode",
            "Chembl smiles",
            "our smiles",
        ],
    )
    smiles_data = pd.DataFrame(smiles_data)
    print(
        f"There are {hetcode_not_found} of {len(combinations_to_investigate)} instances where we skipped a combination because we did not find the hetcode in our smiles, this is {round(hetcode_not_found*100/len(combinations_to_investigate), 2)}%"
    )
    return final_molecule_data, smiles_data, all_hetcodes


def generate_activity_data_for_notfound_entries(notfound_molecules, recalc_all):
    """
    Generates the activity data from the chembl
    :param notfound_molecules:
    :param recalc_all:
    :return:
    """
    if type(notfound_molecules) == pathlib.PosixPath:
        mols_not_found = pd.read_csv(notfound_molecules, index_col=0)
    elif type(notfound_molecules) == pandas.DataFrame:
        mols_not_found = notfound_molecules
    else:
        raise RuntimeError
    number_occurences = []
    all_activities = []
    combinations = {}
    for index, row in mols_not_found.iterrows():
        dataframe = None
        chembl_target = row["Chembl Target"]
        chembl_mol = row["Chembl Molecule"]
        dataframe_path = (
                definitions.unexplainable_mol_path / f"{chembl_target}_{chembl_mol}.csv"
        )
        if recalc_all or not dataframe_path.exists():
            dataframe = query_activities(
                chembl_mol=chembl_mol, chembl_target=chembl_target
            )
            number_occurences.append(dataframe.shape[0])
            dataframe.to_csv(dataframe_path)
        else:
            if dataframe_path.exists():
                dataframe = pd.read_csv(dataframe_path)
                number_occurences.append(dataframe.shape[0])
            else:
                number_occurences.append(0)
        if dataframe is not None:
            for index, row in dataframe.iterrows():
                activity_id = row["activity_id"]
                target_id = row["target_id"]
                molecule_id = row["molecule_id"]
                combinations[activity_id] = [target_id, molecule_id]
                all_activities.append(activity_id)
    mols_not_found["N activities"] = number_occurences
    mols_not_found.to_csv(
        definitions.output_files_path / "smiles_of_molecules_notfound_nactivities.csv"
    )
    return all_activities, combinations, mols_not_found


def split_combinations_by_existing_activities(molecules_not_found_path, recalc=True):
    (
        all_activities,
        combinations,
        mols_not_found,
    ) = generate_activity_data_for_notfound_entries(molecules_not_found_path, recalc)
    # First step done, generate the data of existing activties in the chembl
    # Please note that now in mols_not_found the number of activities in the chembl is annotated, exploit later
    # now calc all activities in the activitydb
    sql_query = f"SELECT * FROM activity WHERE chembl_activity_id IN ({','.join(str(x) for x in all_activities)})"
    activitydb_data = query_activitydb(sql_query)
    relevant_activities = [x for x in activitydb_data["chembl_activity_id"]]
    # Filter the activitycombinations by valid activities
    valid_combinations = {
        x: combinations[x] for x in combinations.keys() if x in relevant_activities
    }
    filtered_combinations = {
        x: combinations[x] for x in combinations.keys() if x not in relevant_activities
    }
    # Split dataframe in ones where the activities are in the activitydb and ones which are not
    # These extra steps are necessary, because the activityIds are not part of the biochemgraph document, only the targets and mols
    with_activities, without_activties = calc_keepers_for_mol_df(
        mols_not_found, valid_combinations
    )
    print(
        f"There are {len(relevant_activities)} activities found for the {len(all_activities)} activities that could possibly be found"
    )
    print(
        f"From the {len(combinations)} combinations found in our data, {len(filtered_combinations)} without activity in our data were found found, which means {round(100*len(valid_combinations)/len(combinations),2)}% remain."
    )
    if len(combinations) != len(filtered_combinations) + len(valid_combinations):
        raise RuntimeError(
            f"We somehow miss combinations when splitting by found activities"
        )
    return (
        valid_combinations,
        filtered_combinations,
        with_activities,
        without_activties,
    )


def calc_keepers_for_mol_df(mols_not_found, valid_combinations):
    with_activities = []
    without_activities = []
    for index, row in mols_not_found.iterrows():
        current_mol = row["Chembl Molecule"]
        current_target = row["Chembl Target"]
        found = False
        for comb in valid_combinations:
            valid_target, valid_mol = valid_combinations[comb]
            if current_mol == valid_mol and current_target == valid_target:
                with_activities.append(index)
                found = True
                break
        if not found:
            without_activities.append(index)
    mols_with_activities = mols_not_found.loc[with_activities]
    mols_without_activities = mols_not_found.loc[without_activities]
    return mols_with_activities, mols_without_activities


def annotate_molecule_df_with_found_chembl_ligands(df_to_annotate):
    found_ligands = []
    did_we_do_it = []
    for index, row in df_to_annotate.iterrows():
        hetcode = row["hetcode"]
        query = f"""
                SELECT *
                FROM pdb_ligand pl
                JOIN pdb_ligand_chembl_ligand cpl ON pl.pdblid = cpl.pdblid
                JOIN chembl_ligand cl ON cpl.clid = cl.clid
                WHERE pl.ligand_id = '{hetcode}';
                """
        chembl_mol = row["Chembl Molecule"]
        dataframe = query_activitydb(query)
        all_found_chembl_ids = [x for x in dataframe["chembl_id"]]
        found_ligands.append(all_found_chembl_ids)
        did_we_do_it.append(chembl_mol in all_found_chembl_ids)

    df_to_annotate["Found ligand"] = did_we_do_it
    df_to_annotate["Chembl ligands"] = found_ligands

    return df_to_annotate


def enrich_and_split_by_chembl_mol_reproducing(df_to_enrich):
    readin_data = get_activitydb_file_processing()
    ligand_data = pd.read_csv(definitions.input_files_path / "parsed_ligand_data.csv")
    valid_mols_not_found = df_to_enrich.merge(
        readin_data, left_on="PDB", right_on="pdb_id", how="left"
    )
    valid_mols_not_found_to_save = valid_mols_not_found.merge(
        ligand_data,
        left_on=["PDB", "hetcode"],
        right_on=["pdb", "hetcodes"],
        how="left",
    )
    valid_mols_not_found_to_save.to_csv(
        definitions.output_files_path
        / "valid_mols_with_found_activities_withchemblligands.csv"
    )

    valid_mols_not_found_to_save_found = valid_mols_not_found_to_save[
        valid_mols_not_found_to_save["Found ligand"] == True
    ]
    valid_mols_not_found_to_save_found.to_csv(
        definitions.output_files_path
        / "valid_mols_with_found_activities_withchemblligands_foundagain.csv"
    )
    valid_mols_not_found_to_save_not_found = valid_mols_not_found_to_save[
        valid_mols_not_found_to_save["Found ligand"] == False
    ]
    valid_mols_not_found_to_save_not_found.to_csv(
        definitions.output_files_path
        / "valid_mols_with_found_activities_withchemblligands_notfoundagain.csv"
    )

    valid_mols_we_did_find = valid_mols_not_found[
        valid_mols_not_found["Found ligand"] == True
    ]
    valid_mols_not_found = valid_mols_not_found[
        valid_mols_not_found["Found ligand"] == False
    ]
    return valid_mols_we_did_find, valid_mols_not_found, valid_mols_not_found_to_save


def get_inchi_from_hetcode(hetcode):
    try:
        url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{hetcode}"
        reply = requests.get(url)
        inchi = reply.json()["rcsb_chem_comp_descriptor"]["in_ch_i"]
    except:
        return None
    return inchi
