import pandas as pd
import json
from src import definitions
from src.utils import db, getter, visualisation


def get_pdbs_in_df(dataframe):
    return len(set(x for x in dataframe["PDB"]))

def get_pdb_target_combinations_in_df(dataframe):
    combinations = set()
    for index, row in dataframe.iterrows():
        pdb = row["PDB"]
        target = row["Target"]
        combinations.add(f"{pdb}-{target}")
    return len(combinations)

recalc = False
pdb_readin_data = getter.get_activitydb_file_processing(recalc=recalc)
ligandfinder_data = pd.read_csv(definitions.input_files_path / "parsed_ligand_data.csv")
with open(
        definitions.output_files_path / "chembiograph.json", "r", encoding="utf-8"
) as file:
    chembiodata_dict = json.load(file)

with open(
        definitions.output_files_path / "activityfinder_target.json", "r", encoding="utf-8"
) as file:
    activityfinder_dict = json.load(file)

with open(
        definitions.output_files_path / "chembiograph_targets.json", "r", encoding="utf-8"
) as file:
    chembio_targets = json.load(file)

chembiodata_set_pdbtarget = set(str(x.upper()) + ":" + str(y) for x in chembio_targets for y in chembio_targets[x])
activityfinder_targetmode_set = set(str(x.upper()) + ":" + str(y) for x in activityfinder_dict for y in activityfinder_dict[x])

visualisation.plot_venn_diagram(chembiodata_set_pdbtarget, activityfinder_targetmode_set, "BioChemGraph", "ActivityFinderTargetMode", definitions.plot_path / "venn_pdb_targetmode_new.png")
visualisation.plot_supervenn_diagram(chembiodata_set_pdbtarget, activityfinder_targetmode_set, "BioChemGraph", "ActivityFinderTargetMode", definitions.plot_path / "venn_super_pdb_targetmode_new.png")

print(len(chembiodata_set_pdbtarget))
print(len(activityfinder_targetmode_set))
print(len(activityfinder_targetmode_set.intersection(chembiodata_set_pdbtarget)))

# Now to the PDBS where we cannot find the correct targets no matter what:
not_found_targets = chembiodata_set_pdbtarget - activityfinder_targetmode_set
not_found_targets = [[x.split(":")[0], x.split(":")[1]] for x in not_found_targets]
notfound_targets = pd.DataFrame(not_found_targets, columns=["PDB", "Target"])



print(notfound_targets.shape)
analysis = notfound_targets.merge(
    pdb_readin_data, left_on="PDB", right_on="pdb_id", how="left"
)
analysis.to_csv(
    definitions.output_files_path
    / "pdb_target_combinations_not_found_targetmode_with_readin.csv"
)
print(analysis.shape)
print(f"We find {analysis.shape[0]} cases where we do not find the correct pdb target combination")
not_allowed = [
    "NMR_SOLUTION",
    "NMR_SOLIDSTATE",
    "NEUTRON_DIFFRACTION",
    "ELECTRON_MICROSCOPY",
    "ELECTRON_CRYST",
]

pdb_target_combinations = {}
analysis_pdb_mirror = analysis[~analysis["comment"].isnull()]
interesting_stuff = analysis_pdb_mirror[~analysis_pdb_mirror["comment"].isin(not_allowed)]
no_readin_issue = analysis_pdb_mirror[analysis_pdb_mirror["comment"] == "success"]
print(f"Of these {analysis.shape[0]}/{analysis.shape[0]/len(chembiodata_set_pdbtarget)} cases, {analysis.shape[0] - analysis_pdb_mirror.shape[0]}/{(analysis.shape[0] - analysis_pdb_mirror.shape[0])/len(chembiodata_set_pdbtarget)} are not in our pdb mirror")
print(f"However, only {no_readin_issue.shape[0]}/{no_readin_issue.shape[0]/analysis.shape[0]} of these cases are without readin issues")
print(analysis["comment"].value_counts())
print(f"There are {interesting_stuff.shape[0] - no_readin_issue.shape[0]} cases where no ligand was found")
for index, row in interesting_stuff.iterrows():
    pdb = row["PDB"]
    target = row["Target"]
    if not pdb in pdb_target_combinations:
        pdb_target_combinations[pdb] = {target: set()}
    else:
        pdb_target_combinations[pdb][target] = set()

bioactivity_data = pd.read_csv(
    definitions.biochemgraph_full, sep="\t"
)
print(bioactivity_data.head())
for index, row in bioactivity_data.iterrows():
    pdb = row["PDB"].upper()
    hetcode = row["LIG"]
    target = row["CHEMBL_ID_TAR"]
    if pdb in pdb_target_combinations:
        if target in pdb_target_combinations[pdb]:
            pdb_target_combinations[pdb][target].add(hetcode)

all_hetcodes = []
for index, row in interesting_stuff.iterrows():
    pdb = row["PDB"]
    target = row["Target"]
    hetcodes = pdb_target_combinations[pdb][target]
    all_hetcodes.append([x for x in hetcodes])

interesting_stuff["hetcodes"] = all_hetcodes

print(interesting_stuff.shape)
print(ligandfinder_data.head())
interesting_stuff = interesting_stuff.merge(
    ligandfinder_data, left_on="PDB", right_on="pdb", how="left"
)
indexes_to_keep = []
for index, row in interesting_stuff.iterrows():
    # print(row["hetcodes_x"], type(row["hetcodes_x"]), row["hetcodes_x"].split(","))
    # hetcodes_found_chembio = row["hetcodes_x"].replace("[","").replace("]","").replace("\'","").split(",")
    hetcodes_found_chembio = row["hetcodes_x"]
    hetcodes_found_ligandfinder = row["hetcodes_y"]
    if hetcodes_found_ligandfinder in hetcodes_found_chembio:
        indexes_to_keep.append(index)

n_pdbs_prior = get_pdb_target_combinations_in_df(interesting_stuff)
print(f"We filter out an additional {interesting_stuff.shape[0]-len(indexes_to_keep)} cases where the hetcodes are different")
interesting_stuff = interesting_stuff.loc[indexes_to_keep]
# all_columns = interesting_stuff.columns.values
interesting_stuff = interesting_stuff.astype(str)
print(f"We filter out an additional {n_pdbs_prior-get_pdb_target_combinations_in_df(interesting_stuff)} cases where the hetcodes are different")
# all_columns = ['PDB', 'Target', 'fid', 'pdb_id', 'read_status', 'comment', 'hetcodes_x', 'pdb', 'skip_reason', 'hetcodes_y', 'nof_residues', 'usmiles', 'SimpleRotB', 'naomi_type', 'formula', 'heavy_atoms', 'mw']
all_columns = [
    "PDB",
    "Target",
    "fid",
    "pdb_id",
    "read_status",
    "comment",
    "hetcodes_x",
    "pdb",
    "skip_reason",
    "hetcodes_y",
    "nof_residues",
    "usmiles",
    "SimpleRotB",
    "naomi_type",
    "formula",
    "heavy_atoms",
    "mw",
]
interesting_stuff = interesting_stuff.drop_duplicates(subset=all_columns, keep="first")
print(interesting_stuff.columns.difference(["name"]).tolist())

all_pdbs = [x for x in interesting_stuff["PDB"]]
all_targets = [x for x in interesting_stuff["Target"]]
activityfinder_target_list = []
found_instead_data = []
for index, pdb in enumerate(all_pdbs):
    chembio_target_description = db.query_chembl(
        f"SELECT * FROM target_dictionary WHERE target_dictionary.chembl_id = '{all_targets[index]}'"
    )
    if chembio_target_description.shape[0] == 0:
        print(f"Skipping pdb {pdb} and target {all_targets[index]}")
        activityfinder_target_list.append([])
        continue
    chembio_organism = chembio_target_description["organism"][0]
    chembio_target_name = chembio_target_description["pref_name"][0]
    chembio_target_type = chembio_target_description["target_type"][0]
    if pdb in activityfinder_dict.keys():
        activityfinder_targets = activityfinder_dict[pdb]
        # Ok somehow I do not understand anymore why I did this part here :thinking:
        for target in activityfinder_targets:
            actfinder_target_description = db.query_chembl(
                f"SELECT * FROM target_dictionary WHERE target_dictionary.chembl_id = '{target}'"
            )
            if actfinder_target_description.shape[0] == 0:
                continue
            actfinder_organism = actfinder_target_description["organism"][0]
            actfinder_target_name = actfinder_target_description["pref_name"][0]
            actfinder_target_type = actfinder_target_description["target_type"][0]
            current_data = [
                pdb,
                all_targets[index],
                chembio_target_name,
                chembio_target_type,
                chembio_organism,
                target,
                actfinder_target_name,
                actfinder_target_type,
                actfinder_organism,
            ]
            found_instead_data.append(current_data)
    else:
        current_data = [
            pdb,
            all_targets[index],
            chembio_target_name,
            chembio_target_type,
            chembio_organism,
            "Nothing found",
            "",
            "",
            "",
        ]
        found_instead_data.append(current_data)
        activityfinder_targets = []
    activityfinder_target_list.append(activityfinder_targets)
print([x for x in interesting_stuff.columns.values])
interesting_stuff.drop(columns=["Unnamed: 0", "fid"], inplace=True)
interesting_stuff.insert(3, "Targets Activityfinder", activityfinder_target_list)

interesting_stuff.to_csv(
    definitions.output_files_path / "interesting_missing_targets_ligand_data.csv"
)

found_instead_data = pd.DataFrame(
    found_instead_data,
    columns=[
        "PDB",
        "Chembio Target",
        "Chembio name",
        "Chembio Type",
        "Chembio organism",
        "Target",
        "Name",
        "Type",
        "Organism",
    ],
)
found_instead_data.to_csv(
    definitions.output_files_path / "chembl_targets_found_instead.csv"
)


data = interesting_stuff
sequence_data = pd.read_csv(
    definitions.output_files_path / "sequence_identity_data_biochem_full.csv"
)

with open(definitions.output_files_path / "chembiograph_targets.json", "r") as f:
    chembio_target_data = json.load(f)
sequence_data_dict = {}
for index, row in sequence_data.iterrows():
    pdb = row["PDB"]
    uniprot = row["Uniprot ID"]
    target = row["ChEMBL Target"]
    identity = row["Best sequence identity"]
    """if not pdb.lower() in chembio_target_data:
        continue
    current_targets = chembio_target_data[pdb.lower()]
    actual_target = None
    for target in current_targets:
        if uniprot in current_targets[target]:
            actual_target = target
    print(uniprot,actual_target,current_targets)"""
    if pdb in sequence_data_dict.keys():
        sequence_data_dict[pdb][target] = identity
    else:
        sequence_data_dict[pdb] = {target: identity}

overlap_data = []
identity_data = []
for index, row in data.iterrows():
    pdb = row["PDB"]
    target = row["Target"]
    current_chembio_targets = [x for x in chembio_targets[pdb.lower()]]
    if pdb in sequence_data_dict.keys() and target in sequence_data_dict[pdb].keys():
        identity_data.append(sequence_data_dict[pdb][target])
    else:
        identity_data.append(None)
    if pdb in activityfinder_dict:
        current_activity_targets = activityfinder_dict[pdb]
        overlap = [x for x in current_chembio_targets if x in current_activity_targets]
        overlap_data.append(
            [len(current_chembio_targets), len(current_activity_targets), len(overlap)]
        )
    else:
        overlap_data.append([len(current_chembio_targets), 0, 0])

overlap_data = pd.DataFrame(
    overlap_data, columns=["N_Chembio", "N_Activity", "Overlap"]
)
data["N_Chembio"] = overlap_data["N_Chembio"]
data["N_Activity"] = overlap_data["N_Activity"]
data["Overlap"] = overlap_data["Overlap"]
data["Identity"] = identity_data

data.to_csv(
    definitions.output_files_path
    / "interesting_missing_targets_ligand_data_overlap.csv"
)

print(f"Finally, there are still {get_pdb_target_combinations_in_df(data)} cases where we do not know what happend")

filtered_data = data[data["comment"] == "success"]

print(f"Finally, there are still {get_pdb_target_combinations_in_df(filtered_data)} cases where we do not know what happend")
filtered_data = filtered_data[filtered_data["Identity"] >= 80.0]
print(f"Only {get_pdb_target_combinations_in_df(filtered_data)} are left after querying for sequence identity")

filtered_data = filtered_data[filtered_data["skip_reason"] == "NotSkipped"]
print(f"Only {get_pdb_target_combinations_in_df(filtered_data)} are left after filtering out ligands that cannot be read")

filtered_data.to_csv(
    definitions.output_files_path / "unexplainable_missing_targets.csv"
)


# _________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________ALAGLYTRPGLYTYRVALTYRASPTYRMETCYSSERPHEALATHRARGARGPHEALAASNHISTYR_______________________________________
# METGLYGLNTHRGLYLYSLYSSERTHRSERVALSERSERLEUARGGLYTHRGLUTHRVALLEUHISASNILEPROTYRMETGLYASPGLUVALLEUASPGLNASPGLYTHRPHEILEGLUGLULEUILELYSASNTYRASPGLYLYSVALHISGLYASPARGGLUCYSGLYPHE____________________________________GLNLYSASNGLNGLUSERASNGLNILEASPLYSGLUSERHISPROPROARGLYSPHELEUTHRGLUGLNGLNLEUPROGLYALALEUPROPROGLUCYSLYSASNTHRGLUMETALAILEASPASNLYSPROCYSGLYPROHISCYSTYRGLNHISLEUGLUGLYALALYSGLUPHEALAALAALALEUTHRALAGLUASNVALGLUGLUSERSERILEILEALAPROALAPROALAGLUASPVALASPTHRLEUGLYGLYGLYGLYSERGLYGLYGLYGLYSERGLYGLYGLYGLYSERALAALAALAASPGLYSERSERASNHISALAGLYTRPGLYTYRVALTYRASPTYRMETCYSSERPHEALATHRARGARGPHEALAASNHISTYRLEULYSTYRVALGLYILEGLUARGGLUMETGLUILEPRO
# MGQTGKKSTSVSSLRGTETVLHNIPYMGDEVLDQDGTFIEELIKNYDGKVHGDRECGF------------QKNQESNQIDKESHPPRKFLTEQQLPGALPPECKNTEMAIDNKPCGPHCYQHLEGAKEFAAALTAENVE___________________________________________________________________________
# MGQTGKKSTSVSSLRGTETVLHNIPYMGDEVLDQDGTFIEELIKNYDGKVHGDRECGFNDDDDDDDGDDPQKDLEDHRDDKESRPPRKFLTEQQLPGALPPECKNTETALDNKPCGPQCYQHLEGAKEFAAALTAERIK___________________________________________________________________________