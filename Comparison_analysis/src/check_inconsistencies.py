import pandas as pd
from pathlib import Path
import definitions

ehmki_data = pd.read_csv("v_presence_quadruplet_enriched.csv", header=None)
ehmki_data.columns = ["PDB", "pdblid", "hetcode", "clid", "chembl_ligand", "tid", "chembl_target", "in_component", "in_variant", "membership"]
print(ehmki_data.head())


quadruple = set()
for chunk in pd.read_csv(definitions.activities_raw, dtype=str, chunksize=1000):
    for index,row in chunk.iterrows():
        pdb = row["pdb.pdb_id"]
        chembl_ligand = row["chembl_ligand.chembl_id"]
        chembl_target = row["target.chembl_id"]
        hetcode = row["pdb_ligand.ligand_id"]
        quadruple.add((pdb, chembl_ligand, chembl_target, hetcode))

new_data = []
for pdb, ligand, target, hetcode in quadruple:
    new_data.append([pdb, ligand,target, hetcode])
new_data = pd.DataFrame(new_data, columns=["PDB", "chembl_ligand", "chembl_target", "hetcode"])

merged_data = ehmki_data.merge(new_data, how ="outer", on=["PDB", "chembl_ligand", "chembl_target", "hetcode"], indicator=True)
print(merged_data.head())
print(merged_data.shape, new_data.shape, ehmki_data.shape)
only_me_data = merged_data[merged_data["_merge"] == "right_only"]
both_data = merged_data[merged_data["_merge"] == "both"]
only_ehmki_data = merged_data[merged_data["_merge"] == "left_only"]
print(only_ehmki_data.shape, only_me_data.shape, both_data.shape)
only_me_data.to_csv("inconsistencies_only_me.csv")
only_ehmki_data.to_csv("inconsistencies_only_ehmki.csv")
both_data.to_csv("inconsistencies_both.csv")