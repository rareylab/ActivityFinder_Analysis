import pandas as pd
import json
import definitions


def add_to_dict(pdb, target_id, mol_id, dictionary):
    if pdb not in dictionary:
        dictionary[pdb] = {target_id: set([mol_id])}
    else:
        if target_id in dictionary[pdb]:
            dictionary[pdb][target_id].add(mol_id)
        else:
            dictionary[pdb][target_id] = set([mol_id])


def clean_dict(dict_to_clean, depth=2):
    for key1 in dict_to_clean:
        if depth == 2:
            for key2 in dict_to_clean[key1]:
                value = list(dict_to_clean[key1][key2])
                dict_to_clean[key1][key2] = value
        elif depth == 1:
            value = list(dict_to_clean[key1])
            dict_to_clean[key1] = value
    return dict_to_clean


def clean_and_save_dict(dict_to_save, output_name, depth=2):
    dict_to_save = clean_dict(dict_to_save, depth=depth)
    with open(definitions.output_files_path / output_name, "w") as outfile:
        json.dump(dict_to_save, outfile)


def extract_target_data(target_data_path):
    target_dict_total = dict()
    target_dict_100 = dict()
    target_dict_95 = dict()
    target_dict_90 = dict()
    for chunk in pd.read_csv(target_data_path, sep=",", chunksize=10000, dtype=str):
        for index, row in chunk.iterrows():
            pdb = row["pdb.pdb_id"]
            tar_chembl_id = row["target.chembl_id"]
            perc_identity = float(row["blast_match.percent_identity"])
            if pdb in target_dict_total.keys():
                target_dict_total[pdb].add(tar_chembl_id)
            else:
                target_dict_total[pdb] = set([tar_chembl_id])
            if perc_identity >= 100:
                if pdb in target_dict_100.keys():
                    target_dict_100[pdb].add(tar_chembl_id)
                else:
                    target_dict_100[pdb] = set([tar_chembl_id])
            if perc_identity >= 95:
                if pdb in target_dict_95.keys():
                    target_dict_95[pdb].add(tar_chembl_id)
                else:
                    target_dict_95[pdb] = set([tar_chembl_id])
            if perc_identity >= 90:
                if pdb in target_dict_90.keys():
                    target_dict_90[pdb].add(tar_chembl_id)
                else:
                    target_dict_90[pdb] = set([tar_chembl_id])
    return target_dict_total, target_dict_100, target_dict_95, target_dict_90


def extract_activity_data(activity_data_path):
    bioactivity_dict_inchi = dict()
    bioactivity_dict_exact = dict()
    bioactivity_dict_good = dict()
    bioactivity_dict_full = dict()
    best_quality_dict = dict()
    for chunk in pd.read_csv(activity_data_path, sep=",", chunksize=10000, dtype=str):
        for index, row in chunk.iterrows():
            pdb = row["pdb.pdb_id"]
            inchikey = row["pdb_ligand.standard_inchi_key"]
            hetcode = row["pdb_ligand.ligand_id"]
            mol_chembl_id = row["chembl_ligand.chembl_id"]
            small_mol_level = row[
                "small_molecule_matching_confidence_level.confidence_level"
            ]
            tar_chembl_id = row["target.chembl_id"]
            perc_identity = row["blast_match.percent_identity"]
            key = f"{pdb}:{tar_chembl_id}:{mol_chembl_id}"
            if pd.isna(pdb) or pd.isna(tar_chembl_id) or pd.isna(mol_chembl_id):
                continue
            if key in best_quality_dict:
                if perc_identity > best_quality_dict[key][0]:
                    best_quality_dict[key][0] = perc_identity
                if small_mol_level > best_quality_dict[key][1]:
                    best_quality_dict[key][1] = small_mol_level
            else:
                best_quality_dict[key] = [perc_identity, small_mol_level]
            if small_mol_level in definitions.inchi_levels:
                add_to_dict(pdb, tar_chembl_id, mol_chembl_id, bioactivity_dict_inchi)
            if small_mol_level in definitions.exact_levels:
                add_to_dict(pdb, tar_chembl_id, mol_chembl_id, bioactivity_dict_exact)
            if small_mol_level in definitions.good_levels:
                add_to_dict(pdb, tar_chembl_id, mol_chembl_id, bioactivity_dict_good)
            add_to_dict(pdb, tar_chembl_id, mol_chembl_id, bioactivity_dict_full)
    return (
        bioactivity_dict_inchi,
        bioactivity_dict_exact,
        bioactivity_dict_good,
        bioactivity_dict_full,
        best_quality_dict,
    )


def extract_activities_total():
    bioactivities_path = definitions.activities_raw
    bioactivities_path_target = definitions.targets_raw
    (
        target_dict_total,
        target_dict_100,
        target_dict_95,
        target_dict_90,
    ) = extract_target_data(bioactivities_path_target)
    (
        bioactivity_dict_inchi,
        bioactivity_dict_exact,
        bioactivity_dict_good,
        bioactivity_dict_full,
        best_quality_dict,
    ) = extract_activity_data(bioactivities_path)
    with open(
            definitions.output_files_path / "activityfinder_quality.json", "w"
    ) as outfile:
        json.dump(best_quality_dict, outfile)
    clean_and_save_dict(bioactivity_dict_inchi, "activityfinder_inchi.json", depth=2)
    clean_and_save_dict(bioactivity_dict_exact, "activityfinder_exact.json", depth=2)
    clean_and_save_dict(bioactivity_dict_good, "activityfinder_good.json", depth=2)
    clean_and_save_dict(bioactivity_dict_full, "activityfinder_full.json", depth=2)
    clean_and_save_dict(target_dict_total, "activityfinder_target.json", depth=1)
    clean_and_save_dict(target_dict_100, "activityfinder_target_100.json", depth=1)
    clean_and_save_dict(
        target_dict_95, "activityfinder_target_95_inclusive.json", depth=1
    )
    clean_and_save_dict(
        target_dict_90, "activityfinder_target_90_inclusive.json", depth=1
    )
    # extract old stuff to make a delta
    (
        bioactivity_dict_inchi,
        bioactivity_dict_exact,
        bioactivity_dict_good,
        bioactivity_dict_full,
        best_quality_dict,
    ) = extract_activity_data(definitions.input_files_path / "parsed_activity_data_raw_state_april25.csv")
    with open(
            definitions.output_files_path / "activityfinder_quality_old.json", "w"
    ) as outfile:
        json.dump(best_quality_dict, outfile)


def main():
    extract_activities_total()


if __name__ == "__main__":
    main()
