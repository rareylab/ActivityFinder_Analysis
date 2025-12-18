from src import definitions
import json
import pandas as pd

activityfinder_target_100_path = (
        definitions.output_files_path / "activityfinder_target_100.json"
)
biochem_target_path = definitions.output_files_path / "chembiograph_targets.json"

with open(activityfinder_target_100_path, "r") as f:
    activityfinder_target_100 = json.load(f)

with open(biochem_target_path, "r") as f:
    biochem_target = json.load(f)


activityfinder_100_only = {
    x: activityfinder_target_100[x]
    for x in activityfinder_target_100
    if not x.lower() in biochem_target
}
print(len(activityfinder_100_only), len(activityfinder_target_100))

interesting_rows = []
n = 0
found_pdb_targets = set()
for chunk in pd.read_csv(
        definitions.input_files_path / "parsed_activity_data_raw_enriched.csv",
    chunksize=1000,
    dtype=str,
):
    for index, row in chunk.iterrows():
        columns = [x for x in chunk.columns.values]
        pdb = row["pdb.pdb_id"]
        percent_identity = row["blast_match.percent_identity"]
        chembl_ligand = row["chembl_ligand.chembl_id"]
        pdb_ligand = row["pdb_ligand.ligand_id"]
        mm_level = row["small_molecule_matching_confidence_level.confidence_level"]
        chembl_target = row["target.chembl_id"]
        pdb_target = pdb + " " + chembl_target
        if (
            pdb in activityfinder_100_only
            and int(mm_level) == 6
            and float(percent_identity) == 100
        ):
            interesting_rows.append([x for x in row])
            if not pdb_target in found_pdb_targets:
                found_pdb_targets.add(pdb_target)
                print(f"Found {n}")
                n += 1

interesting_df = pd.DataFrame(interesting_rows, columns=columns)
interesting_df.to_csv(definitions.automated_path / "Perfect_matches_no_chembio.csv")

interesting_df = interesting_df[interesting_df["Level_pm"] == "Gold"]
interesting_df.to_csv(
    definitions.automated_path / "Perfect_matches_no_chembio_gold.csv"
)
