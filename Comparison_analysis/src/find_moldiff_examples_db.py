import re
from src import definitions
import pandas as pd
from rdkit import Chem
from src.utils import db, visualisation, getter

molecule_difference_path = definitions.general_output_files / "molecule_difference_data"
if not molecule_difference_path.exists():
    molecule_difference_path.mkdir()
molecule_matching_data_path = (
    molecule_difference_path / "molecule_matching_data_activityfinder.csv"
)
molecule_matching_data_path_full = (
    molecule_difference_path / "molecule_matching_data_joined.csv"
)
recalc = False


def get_biochem_combinations(dataframe):
    return set(
        row["CHEMBL_ID_MOL"] + ":" + row["LIG"] for index, row in dataframe.iterrows()
    )


def get_activityfinder_combinations(dataframe):
    return set(
        row["chembl_id"] + ":" + row["ligand_id"] for index, row in dataframe.iterrows()
    )


def print_information_about_hetcode(
    dataframe, hetcode, biochem=True, activityfinder=True
):
    if biochem:
        specfic_df = dataframe[dataframe["LIG"] == hetcode]
    else:
        specfic_df = dataframe[dataframe["ligand_id"] == hetcode]
    if specfic_df.shape[0] > 0:
        print(f"There are {specfic_df.shape[0]} entries of {hetcode} in the data.")
        for index, row in specfic_df.iterrows():
            if biochem:
                print(
                    f"The ligand {hetcode} is linked to {row['CHEMBL_ID_MOL']} by Biochem"
                )
            if activityfinder:
                print(f"The ligand {hetcode} is linked to {row['chembl_id']} by us")
                print(
                    f"This is done with the small molecule matching level of {row['comment']}"
                )
    else:
        print(f"The ligand {hetcode} is not part of the dataframe")


inchi1 = "InChI=1S/C28H31N5O2/c29-17-19-6-8-20(9-7-19)27(34)31-26(14-18-4-2-1-3-5-18)25-16-22(12-13-30-25)21-10-11-24-23(15-21)28(35)33-32-24/h1-5,10-13,15-16,19-20,26H,6-9,14,17,29H2,(H,31,34)(H2,32,33,35)/t19-,20-,26-/m0/s1"
inchi2 = "InChI=1S/C28H31N5O2/c29-17-19-6-8-20(9-7-19)27(34)31-26(14-18-4-2-1-3-5-18)25-16-22(12-13-30-25)21-10-11-24-23(15-21)28(35)33-32-24/h1-5,10-13,15-16,19-20,26H,6-9,14,17,29H2,(H,31,34)(H2,32,33,35)/t19-,20-,26-/m0/s1"


if not molecule_matching_data_path.exists() or recalc:
    sql_statement = """SELECT 
        cl.chembl_id, cl.chembl_unique_smiles AS usmiles_from_chembl, cl.standard_inchi_key AS inchi_key_from_chembl, cl.standard_inchi AS inchi_from_chembl,
        pl.ligand_id, pl.standard_inchi_key as inchi_key_from_pdb, pl.standard_inchi AS inchi_from_pdb, pl.naomi_unique_smiles_chiral AS usmiles_from_pdb,
        smmcl.confidence_level, smmcl.comment
    FROM 
        pdb_ligand AS pl
    INNER JOIN
        pdb_ligand_chembl_ligand AS pcl ON pl.pdblid = pcl.pdblid
    INNER JOIN pdb_ligand_pdb_id AS pdbid ON pl.pdblid = pdbid.pdblid
    INNER JOIN 
    chembl_ligand AS cl ON pcl.clid = cl.clid
    INNER JOIN small_molecule_matching_confidence_level AS smmcl ON pcl.smmclid = smmcl.smmclid;"""
    data = db.query_activitydb(sql_statement)
    data = data.loc[:, ~data.columns.duplicated()]
    data.drop_duplicates(inplace=True)
    hetcode_inchi_dict = {}
    data.to_csv(molecule_matching_data_path)
else:
    data = pd.read_csv(molecule_matching_data_path, index_col=0)
data.drop_duplicates(inplace=True)
# Idee: mache ein ähnliches mapping für moleküle vom BioChemGraph
# Ancshauen was sich unterscheidet
# Was finden wir, aber an einem anderen Level. Warum?
# Was finden die, was wir nicht finden?
# Was finden wir, was die nicht finden? Aufschlüsseln pro level
simple_biochemgraph = definitions.biochemgraph_simple
biochem = pd.read_csv(simple_biochemgraph, sep="\t")
biochem_ligand_only = biochem[["INCHIKEY", "LIG", "CHEMBL_ID_MOL"]]
biochem_ligand_only = biochem_ligand_only.drop_duplicates()
# Here we make venn diagram stuff for Matthias talks
# First, a venn diagram of the hetcodes in their and our approach
combinations_biochem = get_biochem_combinations(biochem_ligand_only)
combinations_activityfinder = get_activityfinder_combinations(data)

hetcodes_biochem = set(biochem_ligand_only["LIG"].to_list())
hetcodes_activityfinder = set(data["ligand_id"].to_list())
visualisation.plot_venn_diagram(
    set1=hetcodes_biochem,
    set2=hetcodes_activityfinder,
    name1="BioChemGraph",
    name2="ActivityFinder",
    output_filename=molecule_difference_path / "venn_hetcodes.png",
)
visualisation.plot_venn_diagram(
    set1=combinations_biochem,
    set2=combinations_activityfinder,
    name1="BioChemGraph",
    name2="ActivityFinder",
    output_filename=molecule_difference_path / "venn_combinations_chembl_hetcode.png",
)
full_data = data.merge(
    biochem_ligand_only,
    left_on="ligand_id",
    right_on="LIG",
    how="outer",
    indicator=True,
)
full_data.to_csv(molecule_matching_data_path_full)

full_data.reset_index(inplace=True)
refound_data = full_data[full_data["_merge"] == "both"]
biochem_only = full_data[full_data["_merge"] == "right_only"]
we_only = full_data[full_data["_merge"] == "left_only"]
refound_combinations_biochem = get_biochem_combinations(refound_data)
refound_combinations_actfindtoll = get_activityfinder_combinations(refound_data)
print(
    f"We find {len(combinations_biochem)} combinations of Biochemgraph and the {len(combinations_activityfinder)} combinations from the Activityfinder with an overlap of {len(combinations_activityfinder.intersection(combinations_biochem))} combinations."
)
print(
    f"When merging with hetcode only, we get {len(refound_combinations_biochem)}/{len(refound_combinations_actfindtoll)} entries where we both find a matching, {len(get_biochem_combinations(biochem_only))} exclusive to the Biochemgraph and {len(get_activityfinder_combinations(we_only))} exclusive to our approach"
)

refound_data.to_csv(molecule_difference_path / "ligand_data_both.csv")
biochem_only.to_csv(molecule_difference_path / "ligand_data_biochem_only.csv")
we_only.to_csv(molecule_difference_path / "ligand_data_we_only.csv")

# Further investigate the data we do not find, maybe a readin problem?
readin_data = getter.get_activitydb_file_processing(recalc=recalc)
valid_chembl_ligands = getter.get_valid_chembl_ligands(recalc=recalc)

chembl_ligands_biochem_only = biochem_only["CHEMBL_ID_MOL"]
print(
    f"There is {len([x for x in chembl_ligands_biochem_only if not x in valid_chembl_ligands])} ligands that are not in our chembl version that the biochemgraph uses"
)
biochem_pdb_ligand = biochem[["PDB", "LIG"]]
biochem_pdb_ligand["PDB"] = [x.upper() for x in biochem_pdb_ligand["PDB"]]
biochem_pdb_ligand_only = biochem_only.merge(biochem_pdb_ligand, how="left", on="LIG")
biochem_only_readin = biochem_pdb_ligand_only.merge(
    readin_data, left_on="PDB", right_on="pdb_id", how="left"
)
biochem_only_readin.to_csv(
    molecule_difference_path / "ligand_data_biochem_only_readin.csv"
)
# Now I want to make another venn diagram for the hetcodes where we both find something for a hetcode
set_activityfinder = set()
set_biochemgraph = set()
indexes_inequal = []
indexes_equal = []
refound_data.reset_index(inplace=True)
for index, row in refound_data.iterrows():
    hetcode = row["ligand_id"]
    chembl_ligand_biochem = row["CHEMBL_ID_MOL"]
    chembl_ligand_us = row["chembl_id"]
    key_activityfinder = hetcode + ":" + chembl_ligand_us
    key_biochem = hetcode + ":" + chembl_ligand_biochem
    set_biochemgraph.add(key_biochem)
    set_activityfinder.add(key_activityfinder)
    if chembl_ligand_biochem != chembl_ligand_us:
        indexes_inequal.append(index)
    else:
        indexes_equal.append(index)
visualisation.plot_venn_diagram(
    set1=set_biochemgraph,
    set2=set_activityfinder,
    name1="BioChemGraph",
    name2="ActivityFinder",
    output_filename=molecule_difference_path / "venn_gemeinsame_hetcodes.png",
)
print(
    f"Going further into the {len(refound_combinations_biochem)}/{len(refound_combinations_actfindtoll)} combinations where we both find a hetcode, there are {len(set_activityfinder - set_biochemgraph)} where only activityfinder finds them and {len(set_biochemgraph - set_activityfinder)} where only biochem finds"
)

"""for index, row in refound_data.iterrows():
    chembl_biochem = row["CHEMBL_ID_MOL"]
    chembl_we = row["chembl_id"]
    if chembl_biochem != chembl_we:
        indexes_inequal.append(index)
    else:
        indexes_equal.append(index)"""

refound_equal = refound_data.loc[indexes_equal]
refound_inequal = refound_data.loc[indexes_inequal]
print(
    f"From the {refound_data.shape[0]} cases where both methods find a matching, in {refound_equal.shape[0]} cases the matchings are identical and in {refound_inequal.shape[0]} they are different"
)
print(
    f"This means that there are {refound_equal.shape[0]-len(combinations_activityfinder.intersection(combinations_biochem))} duplicate cases in these refound cases"
)
refound_equal.to_csv(molecule_difference_path / "refound_equal.csv")
refound_inequal.to_csv(molecule_difference_path / "refound_inequal.csv")
print("For the equal ones we get these matching levels : ")
print(refound_equal["comment"].value_counts())
print(refound_equal["comment"].value_counts(normalize=True))
# unique case analysis
print_information_about_hetcode(refound_equal, "CQD")
print_information_about_hetcode(refound_equal, "5LO")
print_information_about_hetcode(refound_equal, "BRH")
# print_information_about_hetcode(refound_equal, "MYC")
# print_information_about_hetcode(refound_inequal, "MYC", biochem=False)
myc_only = refound_data[refound_data["ligand_id"] == "MYC"]
print(set(x for x in myc_only["chembl_id"]))
print("For the unequal ones we get these matching levels : ")
print(refound_inequal["comment"].value_counts())
hetcodes_found_equal = set([x for x in refound_equal["ligand_id"]])
hetcodes_found_inequal = set([x for x in refound_inequal["ligand_id"]])
hetcodes_found_both = hetcodes_found_equal.intersection(hetcodes_found_inequal)
print(len(hetcodes_found_both), len(hetcodes_found_equal), len(hetcodes_found_inequal))
refound_inequal_true = refound_inequal[
    refound_inequal["ligand_id"].isin(list(hetcodes_found_both))
]
refound_inequal_false = refound_inequal[
    ~refound_inequal["ligand_id"].isin(list(hetcodes_found_both))
]
print(
    f"However, only {refound_inequal_false.shape[0]} of {refound_inequal_false.shape[0] + refound_inequal_true.shape[0]} are different matchings where we do not also reproduce the correct matching"
)
refound_inequal_false.to_csv(
    molecule_difference_path / "refound_data_true_difference.csv"
)
print(
    f"There are {len(set([x for x in refound_inequal_false['ligand_id']]))} where we really find something different."
)

full_data_new = data.merge(
    biochem_ligand_only,
    left_on=["ligand_id", "chembl_id"],
    right_on=["LIG", "CHEMBL_ID_MOL"],
    how="outer",
    indicator=True,
)
activityfinder_only = full_data_new[full_data_new["_merge"] == "left_only"]
combinations_activityfinde_only = get_activityfinder_combinations(activityfinder_only)
print(
    f"Now lets focus on the {len(combinations_activityfinde_only)} combinations exclusive to ActivityFinder"
)
print(
    f"Due to multiple ligand representations, there are {activityfinder_only.shape[0]-len(combinations_activityfinder)} duplicate representations"
)
print(f"So in total there are {activityfinder_only.shape[0]} mappings")
print(activityfinder_only["comment"].value_counts())
print(activityfinder_only["comment"].value_counts(normalize=True))
activityfinder_only.to_csv(molecule_difference_path / "activityfinder_only.csv")
print_information_about_hetcode(activityfinder_only, "LLT", biochem=False)
print_information_about_hetcode(activityfinder_only, "H4B", biochem=False)
print_information_about_hetcode(activityfinder_only, "3BT", biochem=False)
print_information_about_hetcode(activityfinder_only, "MC", biochem=False)

visualisation.visualize_mols_in_df(
    refound_inequal_false, "molecules_truly_different", sort_by="ligand_id"
)
refound_equal_not_identical = refound_equal[refound_equal["confidence_level"] < 5]
visualisation.visualize_mols_in_df(
    refound_equal_not_identical, "refound_equal_not_inchi_match"
)
refound_equal_not_identical.to_csv("refound_data_different_matching_level.csv")

# PDb
# InChI=1S/C22H27NO/c23-18-22(20-11-5-2-6-12-20)15-13-21(14-16-22)24-17-7-10-19-8-3-1-4-9-19/h1-12,21H,13-18,23H2/b10-7+/t21-,22-
# CHEMBL
# InChI=1S/C22H27NO/c23-18-22(20-11-5-2-6-12-20)15-13-21(14-16-22)24-17-7-10-19-8-3-1-4-9-19/h1-12,21H,13-18,23H2/b10-7+/t21-,22-
# our
# InChI=1S/C22H27NO/c23-18-22(20-11-5-2-6-12-20)15-13-21(14-16-22)24-17-7-10-19-8-3-1-4-9-19/h1-12,21H,13-18,23H2/b10-7+/t21-,22-


# unser chembl
# InChI=1S/C22H27NO/c23-18-22(20-11-5-2-6-12-20)15-13-21(14-16-22)24-17-7-10-19-8-3-1-4-9-19/h1-12,21H,13-18,23H2/b10-7+/t21-,22-
# unser pdb
# InChI=1S/C22H27NO/c23-18-22(20-11-5-2-6-12-20)15-13-21(14-16-22)24-17-7-10-19-8-3-1-4-9-19/h1-12,21H,13-18,23H2/b10-7+/t21-,22-

# Analyse the stuff only we find further
# filter for exact matches
we_only_exact = refound_inequal[refound_inequal["confidence_level"] >= 4]
print(we_only_exact.shape)
# Filter out anything that we do not also find the biochem version so it hits heavier and small mols
indexes_to_keep = []
for index, row in we_only_exact.iterrows():
    mol = Chem.MolFromSmiles(row["usmiles_from_pdb"])
    hetcode = row["ligand_id"]
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    print(hetcode, num_heavy_atoms, hetcode in hetcodes_found_both)
    if num_heavy_atoms > 10:
        print("kept")
        indexes_to_keep.append(index)
we_only_exact = we_only_exact.loc[indexes_to_keep]
we_only_exact.to_csv(molecule_difference_path / "overlap_bigmols_different.csv")
we_only_exact = we_only_exact.head(100)
visualisation.visualize_mols_in_df(
    we_only_exact,
    "overlap_bigmols_different",
    skip_their_chembl=False,
    number_of_entries_per_side=10,
    subimgsize=500,
)

we_only_exact = we_only[we_only["confidence_level"] >= 4]
print(we_only_exact.shape)
# Filter out anything that we do not also find the biochem version so it hits heavier and small mols
indexes_to_keep = []
for index, row in we_only_exact.iterrows():
    mol = Chem.MolFromSmiles(row["usmiles_from_pdb"])
    hetcode = row["ligand_id"]
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    print(hetcode, num_heavy_atoms, hetcode in hetcodes_found_both)
    if num_heavy_atoms > 10:
        print("kept")
        indexes_to_keep.append(index)
we_only_exact = we_only_exact.loc[indexes_to_keep]
we_only_exact.to_csv(molecule_difference_path / "we_only_exact_bigmols.csv")
we_only_exact = we_only_exact.head(100)
visualisation.visualize_mols_in_df(
    we_only_exact,
    "weonly_bigmols_different",
    skip_their_chembl=True,
    number_of_entries_per_side=10,
    subimgsize=500,
)

isotope_regex = "\[\d"

we_only_exact = we_only[we_only["confidence_level"] >= 4]
print(we_only_exact.shape)
# Filter out anything that we do not also find the biochem version so it hits heavier and small mols
indexes_to_keep = []
for index, row in we_only_exact.iterrows():
    usmiles_from_pdb = row["usmiles_from_pdb"]
    usmiles_from_chembl = row["usmiles_from_chembl"]
    mol = Chem.MolFromSmiles(usmiles_from_pdb)
    hetcode = row["ligand_id"]
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    print(hetcode, num_heavy_atoms, hetcode in hetcodes_found_both)
    if num_heavy_atoms > 10 and re.search(isotope_regex, usmiles_from_chembl):
        print("kept")
        indexes_to_keep.append(index)
we_only_exact = we_only_exact.loc[indexes_to_keep]
we_only_exact.to_csv(molecule_difference_path / "we_only_exact_bigmols.csv")
we_only_exact = we_only_exact.head(100)
visualisation.visualize_mols_in_df(
    we_only_exact,
    "isotopes",
    skip_their_chembl=True,
    number_of_entries_per_side=10,
    subimgsize=500,
)


we_only_different = we_only[we_only["confidence_level"] < 4]
print(we_only_different.shape)
# Filter out anything that we do not also find the biochem version so it hits heavier and small mols
indexes_to_keep = []
for index, row in we_only_different.iterrows():
    usmiles_from_pdb = row["usmiles_from_pdb"]
    usmiles_from_chembl = row["usmiles_from_chembl"]
    mol = Chem.MolFromSmiles(usmiles_from_pdb)
    hetcode = row["ligand_id"]
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    print(hetcode, num_heavy_atoms, hetcode in hetcodes_found_both)
    if num_heavy_atoms > 10 and re.search(isotope_regex, usmiles_from_chembl):
        print("kept")
        indexes_to_keep.append(index)
we_only_different = we_only_different.loc[indexes_to_keep]
we_only_different.to_csv(molecule_difference_path / "we_only_different_bigmols.csv")
we_only_different = we_only_different.head(100)
visualisation.visualize_mols_in_df(
    we_only_exact,
    "we_only_different",
    skip_their_chembl=True,
    number_of_entries_per_side=10,
    subimgsize=500,
)
