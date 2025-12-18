from pathlib import Path
import os

directory_path = Path(os.path.abspath(__file__)).parent.parent
general_output_files = directory_path / "state_presubmission"
input_files_path = directory_path / "input_files"
output_files_path = general_output_files / "output_files"
output_molecule_path = general_output_files  / "output_molecule_files"
automated_path = general_output_files / "automated_analysis_outputs"
mutation_data_path = general_output_files / "mutation_data"
unexplainable_mol_path = general_output_files / "still_unexplainable_molecules"
molecule_difference_path = general_output_files / "molecule_difference_data"
plot_path = general_output_files / "plots"
extracted_seqres_data = output_files_path / "pdb_seqres_data.json"
all_pdb_sequence_types = ["atomseq", "seqres", "seqresseqadv"]
biochemgraph_full = input_files_path / "bioactivity_report_full_28april.tsv"
biochemgraph_simple = input_files_path / "bioactivity_report_simple_28april.tsv"
activities_raw = input_files_path / "parsed_activity_data_raw_presubmission.csv"
activities_enriched = input_files_path / "parsed_activity_data_raw_enriched_presubmission.csv"
targets_raw = input_files_path / "parsed_target_data_raw_presubmission.csv"
assert targets_raw.exists()
assert activities_raw.exists()
assert activities_enriched.exists()
assert biochemgraph_full.exists()
assert biochemgraph_simple.exists()

ssh_host =   # The name you use in your ssh tunnel configuration
remote_host =   # Hostname of the database server
remote_port =   # Port of the PostgreSQL-Database
db_username =   # activityb username
db_password =   # activitydb password
db_name = "activitydb_test"  # name of the activityDB
db_username2 = "chembl_user"  # name of the chembl user
db_password2 = # chembl user password
db_name2 = "chembl"  # Name of the chembl database

all_levels = [1, 2, 3, 4, 5, 6]
inchi_levels = [5, 6]
exact_levels = [4, 5, 6]
good_levels = [3, 4, 5, 6]

home = False

# Generate all directories if they are not present
general_output_files.mkdir(exist_ok=True)
input_files_path.mkdir(exist_ok=True)
output_files_path.mkdir(exist_ok=True)
output_molecule_path.mkdir(exist_ok=True)
automated_path.mkdir(exist_ok=True)
mutation_data_path.mkdir(exist_ok=True)
unexplainable_mol_path.mkdir(exist_ok=True)
plot_path.mkdir(exist_ok=True)
molecule_difference_path.mkdir(exist_ok=True)
