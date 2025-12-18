from src import definitions
import pandas as pd
import json


def read_mutation_data(mutation_file):
    data = pd.read_csv(mutation_file)
    mutation_dict = {}
    for index, row in data.iterrows():
        pdb_sequence = row["type"]
        assay_id = row["chembl_id"]
        ligand_name = row["naomi_ligand_id"]
        pdb_id = row["pdb_id"]
        sequence_number = row["sequence_number"]
        query_residue = row["blast_query_seq_residue"]
        hit_residue = row["blast_hit_seq_residue"]
        mutation = [sequence_number, query_residue, hit_residue]
        if ligand_name in mutation_dict:
            if assay_id in mutation_dict[ligand_name]:
                if pdb_sequence in mutation_dict[ligand_name][assay_id]:
                    mutation_dict[ligand_name][assay_id][pdb_sequence].append(mutation)
                else:
                    mutation_dict[ligand_name][assay_id][pdb_sequence] = [mutation]
            else:
                mutation_dict[ligand_name][assay_id] = {pdb_sequence: [mutation]}
        else:
            mutation_dict[ligand_name] = {assay_id: {pdb_sequence: [mutation]}}
    return mutation_dict


if __name__ == "__main__":

    complete_mutation_dict = {}

    for file in definitions.mutation_data_path.glob("*.csv"):
        pdb, sequence = file.stem.split("_")
        mutation_data = read_mutation_data(file)
        print(pdb, sequence)
        complete_mutation_dict[pdb] = {sequence: mutation_data}

    with open(
            definitions.output_files_path / "extracted_mutation_data.json", "w"
    ) as outfile:
        json.dump(complete_mutation_dict, outfile)
