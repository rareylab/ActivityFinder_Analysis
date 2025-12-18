import json

import pandas as pd
import requests

import definitions
import sys
from src.utils import db, visualisation


def query_molecules_chembl(chembl_mol_ids):
    molecule_smiles = []
    not_found = 0
    found_mols = []
    for i in range(0, len(chembl_mol_ids), 50):
        current_mols = chembl_mol_ids[i : i + 50]
        cmpd_chembl_ids = ",".join(
            current_mols
        )  # Amend the format of the text string of compounds so that it is suitable for the API call
        limit = 400  # Limit the number of records pulled back for each url call

        # Set up the call to the ChEMBL 'molecule' API
        # Remember that there is a limit to the number of records returned in any one API call (default is 20 records, maximum is 1000 records)
        # So need to iterate over several pages of records to gather all relevant information together!
        url_stem = "https://www.ebi.ac.uk"  # This is the stem of the url
        url_full_string = (
            url_stem
            + "/chembl/api/data/molecule.json?molecule_chembl_id__in={}&limit={}".format(
                cmpd_chembl_ids, limit
            )
        )
        url_full = requests.get(
            url_full_string
        ).json()  # This calls the information back from the API using the 'requests' module, and converts it to json format
        url_molecules = url_full[
            "molecules"
        ]  # This is a list of the results for activities
        while url_full["page_meta"]["next"]:
            url_full = requests.get(url_stem + url_full["page_meta"]["next"]).json()
            url_molecules = (
                url_molecules + url_full["molecules"]
            )  # Add result (as a list) to previous list of results
        for molecule in url_molecules:
            smiles = molecule["molecule_structures"]["canonical_smiles"]
            found_mols.append(molecule["molecule_chembl_id"])
            molecule_smiles.append(smiles)

    print(len(molecule_smiles), len(chembl_mol_ids), len(found_mols))
    not_found += len(chembl_mol_ids) - len(found_mols)
    mol_df = pd.DataFrame({"SMILES": molecule_smiles, "names": found_mols})
    print(f"Did not find {not_found} Molecules")
    return mol_df


pdb_dataframe = db.query_activitydb("SELECT * FROM pdb_file_processing;")
pdb_dataframe.to_csv(
    definitions.output_files_path / "pdbs_processed_activityfinder.csv"
)
# Now we filter to success, which I have to investigate if I want to do that
# pdb_dataframe = pdb_dataframe[pdb_dataframe["comment"] == "success"]
valid_pdbs = set(
    pdb_dataframe["pdb_id"]
)  # All pdbs that we processed with the rolactif
target_df_chembl = db.query_chembl("SELECT * FROM target_dictionary")
# valid_chembl_targets.to_csv(output_files_path / "valid_chembl_targets.csv", index=False)
valid_chembl_targets = [
    x for x in target_df_chembl["chembl_id"]
]  # All chembl targets in our version of chembl

with open(
    definitions.output_files_path / "activityfinder_full.json", "r", encoding="utf-8"
) as file:
    activityfinder_dict = json.load(file)


with open(
    definitions.output_files_path / "chembiograph.json", "r", encoding="utf-8"
) as file:
    chembiodata_dict = json.load(file)


#interesting_from_chembl = pd.read_csv("chembl_info_on_targets_not_found.csv")
#print(interesting_from_chembl)
#ok_targets = [x for x in interesting_from_chembl["chembl_id"]]
print(len(chembiodata_dict))
# Filter out all PDB structures that we did not process in the ActivityDB creation
chembiodata_dict = {
    x.upper(): chembiodata_dict[x]
    for x in chembiodata_dict.keys()
}
print(len(activityfinder_dict))

activityfinder_set_full = set()
chembiodata_set_full = set()
activityfinder_set_pdbtarget = set()
chembiodata_set_pdbtarget = set()
chembiodata_set_pdb = set(x for x in chembiodata_dict.keys())
activityfinder_set_pdb = set(x for x in activityfinder_dict.keys())


for pdb in chembiodata_dict:
    for target in chembiodata_dict[pdb]:
        for ligand in chembiodata_dict[pdb][target]:
            full_entry = pdb + ":" + target + ":" + ligand
            target_entry = pdb + ":" + target
            chembiodata_set_full.add(full_entry)
            chembiodata_set_pdbtarget.add(target_entry)

for pdb in activityfinder_dict:
    for target in activityfinder_dict[pdb]:
        for ligand in activityfinder_dict[pdb][target]:
            full_entry = pdb + ":" + target + ":" + ligand
            target_entry = pdb + ":" + target
            activityfinder_set_full.add(full_entry)
            activityfinder_set_pdbtarget.add(target_entry)


with open(
    definitions.output_files_path / "activityfinder_inchi.json", "r", encoding="utf-8"
) as file:
    activityfinder_dict_inchi = json.load(file)

activityfinder_inchi_set_full = set()
activityfinder_inchi_set_pdbtarget = set()
activityfinder_inchi_set_pdb = set(x for x in activityfinder_dict_inchi.keys())

for pdb in activityfinder_dict_inchi:
    for target in activityfinder_dict_inchi[pdb]:
        for ligand in activityfinder_dict_inchi[pdb][target]:
            full_entry = pdb + ":" + target + ":" + ligand
            target_entry = pdb + ":" + target
            activityfinder_inchi_set_full.add(full_entry)
            activityfinder_inchi_set_pdbtarget.add(target_entry)

with open(
    definitions.output_files_path / "activityfinder_exact.json", "r", encoding="utf-8"
) as file:
    activityfinder_dict_exact = json.load(file)

activityfinder_exact_set_full = set()
activityfinder_exact_set_pdbtarget = set()
activityfinder_exact_set_pdb = set(x for x in activityfinder_dict_exact.keys())

for pdb in activityfinder_dict_exact:
    for target in activityfinder_dict_exact[pdb]:
        for ligand in activityfinder_dict_exact[pdb][target]:
            full_entry = pdb + ":" + target + ":" + ligand
            target_entry = pdb + ":" + target
            activityfinder_exact_set_full.add(full_entry)
            activityfinder_exact_set_pdbtarget.add(target_entry)

with open(
    definitions.output_files_path / "activityfinder_good.json", "r", encoding="utf-8"
) as file:
    activityfinder_dict_good = json.load(file)

activityfinder_good_set_full = set()
activityfinder_good_set_pdbtarget = set()
activityfinder_good_set_pdb = set(x for x in activityfinder_dict_good.keys())

for pdb in activityfinder_dict_good:
    for target in activityfinder_dict_good[pdb]:
        for ligand in activityfinder_dict_good[pdb][target]:
            full_entry = pdb + ":" + target + ":" + ligand
            target_entry = pdb + ":" + target
            activityfinder_good_set_full.add(full_entry)
            activityfinder_good_set_pdbtarget.add(target_entry)

with open(
    definitions.output_files_path / "activityfinder_target.json", "r", encoding="utf-8"
) as file:
    activityfinder_dict_target = json.load(file)

activityfinder_targetmode_set = set()

for pdb in activityfinder_dict_target:
    for target in activityfinder_dict_target[pdb]:
        target_entry = pdb + ":" + target
        activityfinder_targetmode_set.add(target_entry)

# Now to the PDBS where we cannot find the correct targets no matter what:
not_found_targets = chembiodata_set_pdbtarget - activityfinder_targetmode_set
not_found_targets = [[x.split(":")[0], x.split(":")[1]] for x in not_found_targets]
not_found_targets = pd.DataFrame(not_found_targets, columns=["PDB", "Target"])
print(not_found_targets)
print(not_found_targets.head())
not_found_targets.to_csv(
    definitions.output_files_path / "pdb_target_combinations_not_found_targetmode.csv"
)


if not definitions.plot_path.exists():
    definitions.plot_path.mkdir()

visualisation.plot_venn_diagram(
    chembiodata_set_pdbtarget,
    activityfinder_targetmode_set,
    "BioChemGraph",
    "ActivityFinderTargetmode",
    definitions.plot_path / "venn_pdbtarget_targetmode.png",
    shift=True,
)
visualisation.plot_venn_diagram(
    chembiodata_set_pdbtarget,
    activityfinder_set_pdbtarget,
    "BioChemGraph",
    "ActivityFinder",
    definitions.plot_path / "venn_pdbtarget.png",
    "complete",
)
visualisation.plot_venn_diagram(
    chembiodata_set_pdb,
    activityfinder_set_pdb,
    "BioChemGraph",
    "ActivityFinder",
    definitions.plot_path / "venn_pdb.png",
    "complete",
)

visualisation.plot_venn_diagram(
    chembiodata_set_full,
    activityfinder_inchi_set_full,
    "BioChemGraph",
    "ActivityFinder",
    definitions.plot_path / "venn_full_inchi.png",
    "inchi",
)
visualisation.plot_venn_diagram(
    chembiodata_set_pdbtarget,
    activityfinder_inchi_set_pdbtarget,
    "BioChemGraph",
    "ActivityFinder",
    definitions.plot_path / "venn_pdbtarget_inchi.png",
    "inchi",
)
visualisation.plot_venn_diagram(
    chembiodata_set_pdb,
    activityfinder_inchi_set_pdb,
    "BioChemGraph",
    "ActivityFinder",
    definitions.plot_path / "venn_pdb_inchi.png",
    "inchi",
)

visualisation.plot_venn_diagram(
    chembiodata_set_full,
    activityfinder_exact_set_full,
    "BioChemGraph",
    "ActivityFinder",
    definitions.plot_path / "venn_full_exact.png",
    "exact",
)
visualisation.plot_venn_diagram(
    chembiodata_set_pdbtarget,
    activityfinder_exact_set_pdbtarget,
    "BioChemGraph",
    "ActivityFinder",
    definitions.plot_path / "venn_pdbtarget_exact.png",
    "exact",
)
visualisation.plot_venn_diagram(
    chembiodata_set_pdb,
    activityfinder_exact_set_pdb,
    "BioChemGraph",
    "ActivityFinder",
    definitions.plot_path / "venn_pdb_exact.png",
    "exact",
)

visualisation.plot_venn_diagram(
    chembiodata_set_full,
    activityfinder_good_set_full,
    "BioChemGraph",
    "ActivityFinder",
    definitions.plot_path / "venn_full_good.png",
    "good",
)
visualisation.plot_venn_diagram(
    chembiodata_set_pdbtarget,
    activityfinder_good_set_pdbtarget,
    "BioChemGraph",
    "ActivityFinder",
    definitions.plot_path / "venn_pdbtarget_good.png",
    "good",
)
visualisation.plot_venn_diagram(
    chembiodata_set_pdb,
    activityfinder_good_set_pdb,
    "BioChemGraph",
    "ActivityFinder",
    definitions.plot_path / "venn_pdb_good.png",
    "good",
)

if not definitions.output_molecule_path.exists():
    definitions.output_molecule_path.mkdir()

print("Chembio only")
current_subset = chembiodata_set_full - activityfinder_set_full
with open(definitions.output_files_path / "chembio_only_full_all.json", "w") as f:
    json.dump([x for x in current_subset], f)
current_subset = pd.DataFrame(
    [x.split(":") for x in current_subset], columns=["PDB", "Target", "Molecule"]
)
molecules_of_interest = list(set(current_subset["Molecule"]))
molecules_only_chembio = query_molecules_chembl(molecules_of_interest)
molecules_only_chembio.to_csv(
    definitions.output_molecule_path / "molecules_only_chembio.smi",
    sep=" ",
    header=None,
    index=None,
)

print("Chembio only")
current_subset = chembiodata_set_full - activityfinder_exact_set_full
with open(definitions.output_files_path / "chembio_only_full_exact.json", "w") as f:
    json.dump([x for x in current_subset], f)
current_subset = pd.DataFrame(
    [x.split(":") for x in current_subset], columns=["PDB", "Target", "Molecule"]
)
molecules_of_interest = list(set(current_subset["Molecule"]))
molecules_only_chembio = query_molecules_chembl(molecules_of_interest)
molecules_only_chembio.to_csv(
    definitions.output_molecule_path / "molecules_only_chembio_exact.smi",
    sep=" ",
    header=None,
    index=None,
)

print("Acitvity only")
current_subset = activityfinder_exact_set_full - chembiodata_set_full
with open(
    definitions.output_files_path / "activityfinder_only_full_exact.json", "w"
) as f:
    json.dump([x for x in current_subset], f)
current_subset = pd.DataFrame(
    [x.split(":") for x in current_subset], columns=["PDB", "Target", "Molecule"]
)
molecules_of_interest = list(set(current_subset["Molecule"]))
molecules_only_chembio = query_molecules_chembl(molecules_of_interest)
molecules_only_chembio.to_csv(
    definitions.output_molecule_path / "molecules_only_activity_exact.smi",
    sep=" ",
    header=None,
    index=None,
)

print("Both")
current_subset = activityfinder_exact_set_full.union(chembiodata_set_full)
with open(definitions.output_files_path / "both_full_exact.json", "w") as f:
    json.dump([x for x in current_subset], f)
current_subset = pd.DataFrame(
    [x.split(":") for x in current_subset], columns=["PDB", "Target", "Molecule"]
)
molecules_of_interest = list(set(current_subset["Molecule"]))
molecules_only_chembio = query_molecules_chembl(molecules_of_interest)
molecules_only_chembio.to_csv(
    definitions.output_molecule_path / "molecules_both_exact.smi",
    sep=" ",
    header=None,
    index=None,
)

current_subset = chembiodata_set_pdbtarget - activityfinder_set_pdbtarget
with open(definitions.output_files_path / "chembio_only_pdbtarget.json", "w") as f:
    json.dump([x for x in current_subset], f)
