import pandas as pd
import numpy as np
import definitions
from src.utils import visualisation, getter

pdbs_we = set()
pdbs_biochem = set()
ligands_we = set()
ligands_biochem = set()
ligands_we_pdb = set()
unique_ligands_we_pdb = set()
unique_activity_values = set()
ligands_biochem_pdb = set()
targets_we = set()
targets_biochem = set()
combinations_we = set()
combinations_biochem = set()
set_activities = set()

orders = ["PDBs", "ChEMBL Ligands", "ChEMBL Targets", "Combinations"]
sets = [
        [pdbs_we, pdbs_biochem],
        [ligands_we, ligands_biochem],
        [targets_we, targets_biochem],
        [combinations_we, combinations_biochem],
    ]
sets.reverse()
orders.reverse()

# We work with the old calculated data here for ease_of_use
data = pd.read_csv("/work/gutermuth/activityfinderanalysis/src/ownvenn_data.csv")
data['color'] = np.where(data['Label'] == 'ActivityFinder', 'black', 'white')
visualisation.plot_own_venn_diagramm(data, "/work/gutermuth/activityfinderanalysis/state_presubmission/plots/venn_combined.png", "", color=True, order=orders)



print("Activityfinder start")
for chunk in pd.read_csv(
    definitions.activities_raw,
    chunksize=10000,
    dtype=str,
):
    for index, row in chunk.iterrows():
        chembl_ligand = row["chembl_ligand.chembl_id"]
        pdb = row["pdb.pdb_id"]
        target = row["target.chembl_id"]
        hetcode = row["pdb_ligand.ligand_id"]
        inchikey_pdb_ligand = row["pdb_ligand.standard_inchi_key"]
        activity_id = row["activity.chembl_activity_id"]
        combination = pdb.upper() + ":" + target + ":" + chembl_ligand
        unique_activity_values.add(activity_id)
        unique_ligands_we_pdb.add(inchikey_pdb_ligand)
        ligands_we_pdb.add(hetcode)
        pdbs_we.add(pdb)
        targets_we.add(target)
        ligands_we.add(chembl_ligand)
        combinations_we.add(combination)
print("activityfinder done")
print(len(unique_activity_values))
for chunk in pd.read_csv(
    definitions.biochemgraph_full,
    sep="\t",
    chunksize=10000,
    dtype=str,
):
    for index, row in chunk.iterrows():
        pdb = row["PDB"].upper()
        ligand = row["CHEMBL_ID_MOL"]
        target = row["CHEMBL_ID_TAR"]
        hetcode = row["LIG"]
        ligands_biochem.add(ligand)
        ligands_biochem_pdb.add(hetcode)
        targets_biochem.add(target)
        pdbs_biochem.add(pdb)
        combination = pdb.upper() + ":" + target + ":" + ligand
        combinations_biochem.add(combination)


pdb_ligand_data = getter.get_activitydb_pdb_ligand_table(recalc=True)


print(
    f"There are {len(pdbs_we)} PDBS with activity data in our approach and {len(pdbs_biochem)} in the BioChemGraph with an overlap of {len(pdbs_we.intersection(pdbs_biochem))}"
)
print(
    f"There are {len(ligands_we)} CHEMBL Ligands with activity data in our approach and {len(ligands_biochem)} in the BioChemGraph with an overlap of {len(ligands_we.intersection(ligands_biochem))}"
)
print(
    f"There are {len(ligands_we_pdb)} HETCODE ligands with activity data in our approach and {len(ligands_biochem_pdb)} in the BioChemGraph with an overlap of {len(ligands_we_pdb.intersection(ligands_biochem_pdb))}"
)
print(
    f"However, due to the building of multiple molecules per hetcode we get {pdb_ligand_data.shape[0]} different ligands out of the pdb files"
)
print(
    f"There are {len(targets_we)} CHEMBL targets with activity data in our approach and {len(targets_biochem)} in the BioChemGraph with an overlap of {len(targets_we.intersection(targets_biochem))}"
)

visualisation.plot_venn_diagram(
    set1=pdbs_biochem,
    set2=pdbs_we,
    name1="BioChemGraph",
    name2="ActivityFinder",
    output_filename=definitions.plot_path / "venn_pdb_matthias.svg",
    title_addition="PDB IDs",
)
visualisation.plot_venn_diagram(
    set1=ligands_biochem,
    set2=ligands_we,
    name1="BioChemGraph",
    name2="ActivityFinder",
    output_filename=definitions.plot_path / "venn_ligands_matthias.svg",
    title_addition="ChEMBL Ligands",
)
visualisation.plot_venn_diagram(
    set1=ligands_biochem_pdb,
    set2=ligands_we_pdb,
    name1="BioChemGraph",
    name2="ActivityFinder",
    output_filename=definitions.plot_path / "venn_hetcode_matthias.svg",
    title_addition="Hetcode Ligands",
)
visualisation.plot_venn_diagram(
    set1=targets_biochem,
    set2=targets_we,
    name1="BioChemGraph",
    name2="ActivityFinder",
    output_filename=definitions.plot_path / "venn_target_matthias.svg",
    title_addition="ChEMBL Targets",
)
visualisation.plot_venn_diagram(
    set1=combinations_biochem,
    set2=combinations_we,
    name1="BioChemGraph",
    name2="ActivityFinder",
    output_filename=definitions.plot_path / "venn_combination_matthias.svg",
    title_addition="ChEMBL Targets",
)
visualisation.prepare_own_venn_diagramm(
    sets_of_sets=sets,
    labels=["ActivityFinder", "BioChemGraph"],
    output_name=definitions.plot_path / "venn_combined.png",
    title=" ",
    overall_names=orders,
)
