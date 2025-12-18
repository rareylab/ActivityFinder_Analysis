import plotly.graph_objects as go
import pandas as pd
import definitions
import extract_json_from_activityfinder
import extract_json_from_biochemgraph
import json
from src.utils import general, getter, visualisation

def convert_sankey_to_sunburst(dict_to_convert):

    values = dict_to_convert["values"]
    labels = dict_to_convert["labels"]
    sources = dict_to_convert["sources"]
    targets = dict_to_convert["targets"]
    s_labels = ["All combinations"]
    s_parents = [""]
    s_ids = ["0"]
    s_values = [values[0] + values[1] + values[2]]

    i = 0
    for s, t, v in zip(sources, targets, values):
        i += 1
        # if i > 16:
        #    continue
        print(s, t, v, labels[t])
        s_ids.append(str(t))
        s_parents.append(str(s))
        s_values.append(v)
        s_labels.append(labels[t])

    return s_labels, s_parents, s_values, s_ids


# If this is True, everything is rerun starting from the input data
rerun_all = False

# We immediately start to do construct the data for the sankey diagram
global current_highest_sankey_target
global sankey_dict
sankey_dict = {
    "labels": ["All combinations"],
    "sources": [],
    "targets": [],
    "values": [],
}
current_highest_sankey_target = 1


def add_entry_to_sankey_dict(label, source, target, value):
    sankey_dict["labels"].append(label)
    sankey_dict["sources"].append(source)
    sankey_dict["targets"].append(target)
    sankey_dict["values"].append(value)


def split_by_novelty(combinations_to_split, combinations_to_compare):
    # Now we really looked deep into the problems with the ActivityFinder.
    # How bout we now get into the positives?
    # IDEA: We have combinations of PDB (A) ChEMBL Target (B) and ChEMBL Molecule (C)
    # Sort the only activityfinder stuff into novel A / B / C / AB / AC / BC / ABC

    all_pdbs_compare = {pdb for pdb in combinations_to_compare}
    all_targets_compare = {
        target
        for pdb in combinations_to_compare
        for target in combinations_to_compare[pdb]
    }
    all_mols_compare = {
        molecule
        for pdb in combinations_to_compare
        for target in combinations_to_compare[pdb]
        for molecule in combinations_to_compare[pdb][target]
    }
    # Initialize counters
    counts = {
        "different_comb": 0,
        "a": 0,
        "b": 0,
        "c": 0,
        "ab": 0,
        "ac": 0,
        "bc": 0,
        "abc": 0,
    }

    total = 0

    for key in combinations_to_split:
        pdb, target, molecule = key.split(":")
        total += 1

        in_pdb = pdb in all_pdbs_compare
        in_target = target in all_targets_compare
        in_molecule = molecule in all_mols_compare

        if in_pdb and in_target and in_molecule:
            counts["different_comb"] += 1
        elif in_pdb and in_target:
            counts["c"] += 1
        elif in_pdb and in_molecule:
            counts["b"] += 1
        elif in_target and in_molecule:
            counts["a"] += 1
        elif in_pdb:
            counts["bc"] += 1
        elif in_target:
            counts["ac"] += 1
        elif in_molecule:
            counts["ab"] += 1
        else:
            counts["abc"] += 1

    # Check the total counts
    if total != sum(counts.values()) or total != len(combinations_to_split):
        raise RuntimeError("Something went wrong counting here")
    return counts


def split_and_add_by_novelty(set_to_split, set_to_compare, start_point):

    global current_highest_sankey_target
    counts = split_by_novelty(set_to_split, set_to_compare)

    current_highest_sankey_target += 1
    add_entry_to_sankey_dict(
        "Different combination",
        start_point,
        current_highest_sankey_target,
        counts["different_comb"],
    )
    current_highest_sankey_target += 1
    add_entry_to_sankey_dict(
        "Novel PDB", start_point, current_highest_sankey_target, counts["a"]
    )
    current_highest_sankey_target += 1
    add_entry_to_sankey_dict(
        "Novel Target", start_point, current_highest_sankey_target, counts["b"]
    )
    current_highest_sankey_target += 1
    add_entry_to_sankey_dict(
        "Novel Molecule", start_point, current_highest_sankey_target, counts["c"]
    )
    current_highest_sankey_target += 1
    add_entry_to_sankey_dict(
        "Novel PDB and Target", start_point, current_highest_sankey_target, counts["ab"]
    )
    current_highest_sankey_target += 1
    add_entry_to_sankey_dict(
        "Novel PDB and Molecule",
        start_point,
        current_highest_sankey_target,
        counts["ac"],
    )
    current_highest_sankey_target += 1
    add_entry_to_sankey_dict(
        "Novel Target and Molecule",
        start_point,
        current_highest_sankey_target,
        counts["bc"],
    )
    current_highest_sankey_target += 1
    add_entry_to_sankey_dict(
        "Novel everything", start_point, current_highest_sankey_target, counts["abc"]
    )


def split_and_add_by_quality_levels(
    set_to_split, quality_dict, start_point, calc_further=False
):
    # MM level : Perfect, achiral, rest
    # Sequence level : 100%, 95%, rest
    # Initialize counters
    global current_highest_sankey_target
    counts = {
        "perfect_100": [0, []],
        "perfect_95": [0, []],
        "perfect_rest": [0, []],
        "achiral_100": [0, []],
        "achiral_95": [0, []],
        "achiral_rest": [0, []],
        "rest_100": [0, []],
        "rest_95": [0, []],
        "rest_rest": [0, []],
    }

    total = 0

    for key in set_to_split:
        total += 1
        sequence_quality, molecule_quality = quality_dict[key]
        molecule_quality = int(molecule_quality)
        sequence_quality = float(sequence_quality)
        if molecule_quality >= 4:
            # Perfect MM
            if sequence_quality >= 100:
                counts["perfect_100"][0] += 1
                counts["perfect_100"][1].append(key)
            elif sequence_quality >= 95:
                counts["perfect_95"][0] += 1
                counts["perfect_95"][1].append(key)
            else:
                counts["perfect_rest"][0] += 1
                counts["perfect_rest"][1].append(key)
        elif molecule_quality == 3:
            # achiral
            if sequence_quality >= 100:
                counts["achiral_100"][0] += 1
                counts["achiral_100"][1].append(key)
            elif sequence_quality >= 95:
                counts["achiral_95"][0] += 1
                counts["achiral_95"][1].append(key)
            else:
                counts["achiral_rest"][0] += 1
                counts["achiral_rest"][1].append(key)
        else:
            # rest
            if sequence_quality >= 100:
                counts["rest_100"][0] += 1
                counts["rest_100"][1].append(key)
            elif sequence_quality >= 95:
                counts["rest_95"][0] += 1
                counts["rest_95"][1].append(key)
            else:
                counts["rest_rest"][0] += 1
                counts["rest_rest"][1].append(key)

    for key in counts.keys():
        current_highest_sankey_target += 1
        add_entry_to_sankey_dict(
            key, start_point, current_highest_sankey_target, counts[key][0]
        )
        if calc_further:
            split_and_add_by_novelty(
                counts[key][1], chembiodata_dict, current_highest_sankey_target
            )
    return counts


# First step is extracting the keys from the input data
# First the BioChemgraph Data
bioactivities_path = definitions.biochemgraph_full
chembio_dict_path = definitions.output_files_path / "chembiograph.json"
chembio_target_dict = definitions.output_files_path / "chembiograph_targets.json"
if (
    not chembio_dict_path.exists()
    or not chembio_target_dict.exists()
    or rerun_all is True
):
    (
        bioactivity_dict,
        target_uniprot_dict,
    ) = extract_json_from_biochemgraph.extract_bioactivity_report(bioactivities_path)
    extract_json_from_biochemgraph.process_and_save_dicts(
        bioactivity_dict, target_uniprot_dict
    )

with open(
        definitions.output_files_path / "chembiograph.json", "r", encoding="utf-8"
) as file:
    chembiodata_dict = json.load(file)

# We have to do this to make sure that the PDB codes are written identically
chembiodata_dict = {x.upper(): chembiodata_dict[x] for x in chembiodata_dict.keys()}

# Now the ActivityFinder Data
activityfinder_dict_path = definitions.output_files_path / "activityfinder_full.json"
activityfinder_target_dict = (
        definitions.output_files_path / "activityfinder_target.json"
)
activityfinder_quality_path = (
        definitions.output_files_path / "activityfinder_quality.json"
)
bioactivities_path = definitions.activities_raw
bioactivities_path_target = definitions.targets_raw
if (
    not activityfinder_dict_path.exists()
    or not activityfinder_target_dict.exists()
    or not activityfinder_quality_path.exists()
    or rerun_all is True
):
    extract_json_from_activityfinder.extract_activities_total()

with open(activityfinder_dict_path, "r", encoding="utf-8") as file:
    activityfinder_dict = json.load(file)
with open(activityfinder_quality_path, "r", encoding="utf-8") as file:
    activityfinder_quality_dict = json.load(file)

print(
    len(activityfinder_dict.keys()),
    len(chembiodata_dict.keys()),
    len(activityfinder_quality_dict.keys()),
)
# Now we want to get all the data we will need down the line, like which pdbs we processed etc.
# First, the pdb file processing dataframe from the ActivityDB to filter for valid PDB files
"""file_processing_df = query_activitydb("SELECT * FROM pdb_file_processing;")
file_processing_df.to_csv(
    definitions.output_files_path / "pdbs_processed_activityfinder.csv"
)"""
file_processing_df = getter.get_activitydb_file_processing(rerun_all)
processed_pdbs = getter.get_processed_pdbs(
    rerun_all
)  # All pdbs that we processed with the rolactif
valid_pdbs = getter.get_valid_pdbs(
    rerun_all
)  # ALL pdbs that we correctly processed with the rolactif

# All chembl targets in our version of chembl
valid_chembl_targets = getter.get_valid_chembl_targets(rerun_all)
print(
    f"There are {len(valid_chembl_targets)} chembl targets in our chembl version and {len(processed_pdbs)} pdbs out of which {len(valid_pdbs)} {round(len(valid_pdbs) / len(processed_pdbs),4)*100}% were correctly processed"
)
# Now we start comparing stuff

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

chembio_only = chembiodata_set_full - activityfinder_set_full
activity_only = activityfinder_set_full - chembiodata_set_full
union = activityfinder_set_full.union(chembiodata_set_full)
intersection = activityfinder_set_full.intersection(chembiodata_set_full)

visualisation.plot_venn_diagram(
    chembiodata_set_full,
    activityfinder_set_full,
    "BioChemGraph",
    "ActivityFinder",
    definitions.automated_path / "venn_full.svg",
    "complete",
)

visualisation.plot_supervenn_diagram(
    chembiodata_set_full,
    activityfinder_set_full,
    "BioChemGraph",
    "ActivityFinder",
    definitions.automated_path / "supervenn_full.png",
    "complete",
)
n_activityfinder_only = len(activity_only)
n_chembioonly = len(chembio_only)
n_both = len(union)
n_intersection = len(intersection)

print(
    f"There are a total of {n_both} different complexes that we can find ActivityData for in the Chembl"
)
print(
    f"Of these cases, {n_intersection} are found by both approaches, the ActivityFinder finds {n_activityfinder_only} complexes that BioChemGraph does not find and BioChemGraph finds {n_chembioonly} that the ActivityFinder does not find"
)
print(
    f"For now, we want to focus on the {n_chembioonly} cases that we cannot currently reproduce"
)

# Now we can start to populate the sankey data
add_entry_to_sankey_dict("ActivityFinder", 0, 1, n_activityfinder_only)
add_entry_to_sankey_dict("Both", 0, 2, n_intersection)
add_entry_to_sankey_dict("BioChemGraph", 0, 3, n_chembioonly)

# First reason : We do not correctly read in the PDB files
chembio_only_readable = set([x for x in chembio_only if x.split(":")[0] in valid_pdbs])
complexes_filtered_by_readability = chembio_only - chembio_only_readable
n_filtered_readin = n_chembioonly - len(chembio_only_readable)
if n_filtered_readin != len(chembio_only) - len(chembio_only_readable):
    raise RuntimeError("Readin seems malfunictioning")
# Add this as a reason to the sankey dict
add_entry_to_sankey_dict("PDB readin problem", 3, 4, n_filtered_readin)
print(
    f"If we now remove all PDBs that we did not correctly processed, we find the explanation for {n_filtered_readin} complexes, which is {round((n_filtered_readin / n_chembioonly)*100, 2)} % "
)
complexes_filtered = pd.DataFrame(
    [x.split(":") for x in complexes_filtered_by_readability],
    columns=["PDB", "Target", "Molecule"],
)
complexes_filtered = complexes_filtered.merge(
    file_processing_df, left_on="PDB", right_on="pdb_id", how="left"
)
complexes_filtered.to_csv(definitions.output_files_path / "complexes_without_working_readin.csv")
print("There are multiple possible reasons why A PDB cannot be read in.")
specific_readin_dict = {}
PDB_readin_details = complexes_filtered["comment"].value_counts()
complete_count = 0
readin_with_reason = 0
cutoff_showing = 5
current_highest_sankey_target = 4
for name, count in PDB_readin_details.items():
    complete_count += count
    if count < cutoff_showing:
        readin_with_reason += count
    else:
        print(
            f"Due to the problem {name} we cannot reproduce {count} cases, which is {round(count*100 / n_filtered_readin, 2)}% of the filtered cases"
        )
        specific_readin_dict[name] = count
        current_highest_sankey_target += 1
        add_entry_to_sankey_dict(name, 4, current_highest_sankey_target, count)
print(
    f"There are {readin_with_reason} further cases that are filtered with reasons occuring less than {cutoff_showing} times"
)
current_highest_sankey_target += 1
add_entry_to_sankey_dict(
    "Other readin reasons", 4, current_highest_sankey_target, readin_with_reason
)
print(
    f"There are {n_filtered_readin - complete_count} cases ({round((n_filtered_readin - complete_count)*100 / n_filtered_readin, 2)}%) without any reasons, which means that they were not in the used PDB mirror"
)
print(f"There are {complexes_filtered[complexes_filtered["pdb_id"].isna()].shape}")
specific_readin_dict["Other reasons"] = readin_with_reason
specific_readin_dict["PDB mirror"] = n_filtered_readin - complete_count
current_highest_sankey_target += 1
not_in_mirror = n_filtered_readin - complete_count
add_entry_to_sankey_dict(
    "Not in PDB Mirror", 4, current_highest_sankey_target, not_in_mirror
)
if complete_count + not_in_mirror != len(complexes_filtered):
    complexes_filtered.to_csv("debug.csv")
    raise RuntimeError(complete_count, n_filtered_readin, len(complexes_filtered))

# Calc valid chembl ligands
valid_ligands = getter.get_valid_chembl_ligands(rerun_all)
# Now we go into the intricacies of versioning differences
chembio_only_target_found = set(
    [x for x in chembio_only_readable if x.split(":")[1] in valid_chembl_targets]
)
n_filtered_targetfound = len(chembio_only_readable) - len(chembio_only_target_found)
print(
    "Now there are some versioning differences that we should take into account. It might be the case that either the targets, molecules or activities do not exist in our version of ChEMBL."
)
print(
    f"If we now remove all PDBs that we do not have the chembl target, we find the explanation for {n_filtered_targetfound} complexes, which is {round((n_filtered_targetfound / len(chembio_only_readable))*100, 2)} % "
)

chembio_only_mol_found = set(
    [x for x in chembio_only_target_found if x.split(":")[2] in valid_ligands]
)
n_filtered_molfound = len(chembio_only_target_found) - len(chembio_only_mol_found)
print(
    f"If we now remove all PDBs that we do not have the chembl ligand, we find the explanation for {n_filtered_molfound} complexes, which is {round((n_filtered_molfound / len(chembio_only_target_found))*100, 2)} % "
)
if len(chembio_only_mol_found) + n_filtered_molfound + n_filtered_targetfound != len(
    chembio_only_readable
):
    raise RuntimeError("Initial Compatibility seems losing some stuff")
(
    combinations_with_target,
    combinations_without_target,
) = general.split_combinations_by_target_found(chembio_only_mol_found, rerun_all)
#combinations_without_target.to_csv("test_data_without_target.csv")
if len(combinations_with_target) + len(combinations_without_target) != len(
    chembio_only_mol_found
):
    print(
        len(combinations_with_target),
        len(combinations_without_target),
        len(chembio_only_mol_found),
    )
    raise RuntimeError("This split function does not work correctly")
print(
    f"In {len(combinations_without_target)} cases we cannot correctly link the right target"
)


print(f"Combinations with target {len(combinations_with_target)}")
print("Now we need to get the hetcodes and the respective smiles from us")
(
    final_molecule_data,
    smiles_data,
    hetcodes_of_combinations,
) = general.generate_smiles_for_combinations(
    combinations_with_target, recalc_all=False
)
if final_molecule_data.shape[0] != len(combinations_with_target):
    print(final_molecule_data.shape[0], len(combinations_with_target))
    raise RuntimeError("Hier ist ein fehler officer")
(
    combinations_working_adb,
    combinations_not_working_adb,
    with_activities_adb,
    without_activities_adb,
) = general.split_combinations_by_existing_activities(
    final_molecule_data, recalc=False
)
if (
    with_activities_adb.shape[0] + without_activities_adb.shape[0]
    != final_molecule_data.shape[0]
):
    print(
        with_activities_adb.shape[0],
        without_activities_adb.shape,
        final_molecule_data.shape,
    )
    raise RuntimeError("Nein hier ist der Fehler officer")
n_no_activity_in_chembl = without_activities_adb[
    without_activities_adb["N activities"] == 0
].shape[0]
n_no_activity_in_activitydb = without_activities_adb.shape[0] - n_no_activity_in_chembl
without_activities_adb.to_csv(definitions.output_files_path / "data_without_activities_adb.csv")
with_activities_adb.to_csv(definitions.output_files_path / "data_wit_activities_adb.csv")

print(
    f"There are {n_no_activity_in_chembl} combinations where no activity value is in the chembl"
)
print(
    f"There are {n_no_activity_in_activitydb} combinations that are filtered, because no activity ID that is present for BioChemGraph is in the ActivityDB"
)
if (
    n_no_activity_in_chembl + n_no_activity_in_activitydb
    != without_activities_adb.shape[0]
):
    raise RuntimeError("huiarenuiratde")

print(without_activities_adb.columns.values)
print(without_activities_adb.head())
all_chembl_ligands_in_activityb = getter.get_chembl_ligands_in_activitydb()
all_chembl_targets_in_activitydb = getter.get_chembl_targets_in_activitydb()
n_no_activity = 0
n_no_activity_no_ligand = 0
n_no_activity_no_target = 0
n_no_activity_no_ligandandtarget = 0
for index, row in without_activities_adb.iterrows():
    ligand_id = row["Chembl Molecule"]
    target = row["Chembl Target"]
    n_activities = row["N activities"]
    ligand_found = ligand_id in all_chembl_ligands_in_activityb
    target_found = target in all_chembl_targets_in_activitydb
    if n_activities == 0:
        continue
    if target_found and ligand_found:
        n_no_activity += 1
    elif target_found:
        n_no_activity_no_ligand += 1
    elif ligand_found:
        n_no_activity_no_target += 1
    else:
        n_no_activity_no_ligandandtarget += 1

print(f"Of these cases, {n_no_activity} are just without an activity, in {n_no_activity_no_target} the target is not there, in {n_no_activity_no_ligand} the ligand is not there and in {n_no_activity_no_ligandandtarget} both ligand and target are missing")
current_highest_sankey_target += 1
add_entry_to_sankey_dict(
    "Target was not reproduced",
    3,
    current_highest_sankey_target,
    len(combinations_without_target),
)
# Hier könnte ich nochmal nachschauen nach der Sequenzidentität und das weiter aufsplitten
sequence_info_data = getter.get_sequence_identity_data_biochem()
all_targets_not_reproduces = [x[1] for x in combinations_without_target]
keepers = []
combinations_sequence = {}
for index, row in sequence_info_data.iterrows():
    pdb = row["PDB"]
    target = row["ChEMBL Target"]
    for combination in combinations_without_target:
        key = f"{combination[0]}:{combination[1]}:{combination[2]}"
        if (
            pdb == combination[0]
            and target == combination[1]
            and key not in combinations_sequence
        ):
            sequence_identity = row["Best sequence identity"]
            combinations_sequence[key] = sequence_identity
            break

print(len(combinations_without_target), len(combinations_sequence))
low_seq = [x for x in combinations_sequence if combinations_sequence[x] < 80]
not_low_seq = [x for x in combinations_sequence if not combinations_sequence[x] < 80]
print(len(low_seq), len(not_low_seq), len(combinations_without_target))
target_not_reproduced_low_identity = definitions.output_files_path / "target_not_reproduced_lowseqid.json"
target_not_reproduced_unknown = definitions.output_files_path / "target_not_reproduced_unknown.json"
with open(target_not_reproduced_low_identity, "w") as f:
    json.dump(low_seq, f)
with open(target_not_reproduced_unknown, "w") as f:
    json.dump(not_low_seq, f)
print(len(low_seq), len(combinations_without_target) - len(low_seq))
current_highest_sankey_target += 1
add_entry_to_sankey_dict(
    "Sequence identity below 80%",
    current_highest_sankey_target - 1,
    current_highest_sankey_target,
    len(low_seq),
)
current_highest_sankey_target += 1
add_entry_to_sankey_dict(
    "Rest",
    current_highest_sankey_target - 2,
    current_highest_sankey_target,
    len(combinations_without_target) - len(low_seq),
)

n_compatibility_issue = (
    n_no_activity_in_chembl
    + n_no_activity_in_activitydb
    + n_filtered_molfound
    + n_filtered_targetfound
)
current_highest_sankey_target += 1
add_entry_to_sankey_dict(
    "Compatability issues", 3, current_highest_sankey_target, n_compatibility_issue,
)
compatibility_issue_level = current_highest_sankey_target
current_highest_sankey_target += 1
add_entry_to_sankey_dict(
    "Activity not in ChEMBL",
    compatibility_issue_level,
    current_highest_sankey_target,
    n_no_activity_in_chembl,
)
current_highest_sankey_target += 1
add_entry_to_sankey_dict(
    "Activity not in ActivityDB",
    compatibility_issue_level,
    current_highest_sankey_target,
    n_no_activity_in_activitydb,
)
# Hier muss ich noch nachschauen, ob die Aktivität es in unserer CHEMBL version gibt
further_save = current_highest_sankey_target
current_highest_sankey_target += 1
add_entry_to_sankey_dict(
    "Activity with ligand and target in ActivityDB",
    further_save,
    current_highest_sankey_target,
    n_no_activity,
)
current_highest_sankey_target += 1
add_entry_to_sankey_dict(
    "Activity without ligand in ActivityDB",
    further_save,
    current_highest_sankey_target,
    n_no_activity_no_ligand,
)
current_highest_sankey_target += 1
add_entry_to_sankey_dict(
    "Activity without target in ActivityDB",
    further_save,
    current_highest_sankey_target,
    n_no_activity_no_target,
)
current_highest_sankey_target += 1
add_entry_to_sankey_dict(
    "Activity without target and molecule in ActivityDB",
    further_save,
    current_highest_sankey_target,
    n_no_activity_no_ligandandtarget,
)
current_highest_sankey_target += 1
add_entry_to_sankey_dict(
    "Chembl Target not in ChEMBL",
    compatibility_issue_level,
    current_highest_sankey_target,
    n_filtered_targetfound,
)
current_highest_sankey_target += 1
add_entry_to_sankey_dict(
    "Chembl Molecule not in ChEMBL",
    compatibility_issue_level,
    current_highest_sankey_target,
    n_filtered_molfound,
)


# Now we get into the nitty gritty of reproducing the hetcode - chembl ligand partners
annotated_df_path = (
        definitions.output_files_path / "valid_mols_with_activities_with_liganddata.csv"
)
if not annotated_df_path.exists():
    annotated_df = general.annotate_molecule_df_with_found_chembl_ligands(
        with_activities_adb
    )
    annotated_df.to_csv(annotated_df_path)
else:
    annotated_df = pd.read_csv(annotated_df_path)
if annotated_df.shape[0] != with_activities_adb.shape[0]:
    raise RuntimeError("DF annotation seems to change the number of rows...")

(
    valid_mols_we_did_find,
    valid_mols_not_found,
    valid_mols_not_found_to_save,
) = general.enrich_and_split_by_chembl_mol_reproducing(annotated_df)

print(
    f"From the remaining {annotated_df.shape[0]} cases, in {valid_mols_we_did_find.shape[0]} cases we correctly reproduce the chembl ligand and the issue is something with ligand recognition"
)
print(
    f"and in the other {valid_mols_not_found.shape[0]} cases, we do not correlate the hetcode to the same ChEMBL ligand as BioChemGraph."
)

n_ligand_connection = valid_mols_not_found.shape[0]
n_ligand_recognition = valid_mols_we_did_find.shape[0]
n_ligand_not_reproduced = n_ligand_connection + n_ligand_recognition
if annotated_df.shape[0] != n_ligand_not_reproduced:
    raise RuntimeError("Enrich and split my chembl mol seems to forget stuff")

current_highest_sankey_target += 1
add_entry_to_sankey_dict(
    "Ligand was not reproduced",
    3,
    current_highest_sankey_target,
    n_ligand_not_reproduced,
)
current_highest_sankey_target += 1
add_entry_to_sankey_dict(
    "Ligand recognition problem",
    current_highest_sankey_target - 1,
    current_highest_sankey_target,
    n_ligand_recognition,
)
current_highest_sankey_target += 1
add_entry_to_sankey_dict(
    "Hetcode - Chembl connection problem",
    current_highest_sankey_target - 2,
    current_highest_sankey_target,
    n_ligand_connection,
)


total = (
    n_filtered_readin
    + len(combinations_without_target)
    + n_ligand_not_reproduced
    + n_compatibility_issue
)
if total != len(chembio_only):
    raise RuntimeError("Losing some chembio only stuff somewhere")
# split_and_add_by_novelty(activity_only, chembiodata_dict, 1)
# split_and_add_by_novelty(union, chembiodata_dict, 2)
counts_activity = split_and_add_by_quality_levels(
    activity_only, activityfinder_quality_dict, 1, True
)
counts_intersection = split_and_add_by_quality_levels(
    intersection, activityfinder_quality_dict, 2, calc_further=False
)
# split_and_add_by_novelty(counts_activity["perfect_100"][1])

fig = go.Figure(
    data=[
        go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=sankey_dict["labels"],
                x=[0, 0.3, 0.3, 0.3],
                y=[0.5, 0.3, 0.8, 1],
            ),
            link=dict(
                source=sankey_dict["sources"],
                target=sankey_dict["targets"],
                value=sankey_dict["values"],
            ),
        )
    ]
)

fig.update_layout(title_text="BioChemGraph vs ActivityFinder", font_size=10)
fig.write_html(definitions.automated_path / "sankey_diagram.html")
fig.write_image(definitions.automated_path / "sankey_diagram.png")
with open(definitions.output_files_path / "sankey_data.json", "w") as f:
    json.dump(sankey_dict, f)


(
    sunburst_labels,
    sunburst_parents,
    sunburst_values,
    sunburst_ids,
) = convert_sankey_to_sunburst(sankey_dict)

# Erstellung des Sunburst-Diagramms
fig = go.Figure(
    go.Sunburst(
        ids=sunburst_ids,
        labels=sunburst_labels,
        parents=sunburst_parents,
        values=sunburst_values,
        branchvalues="total",
    )
)
fig.update_layout(margin=dict(t=0, l=0, r=0, b=0))
fig.write_html(definitions.automated_path / "sunburst_new.html")
