import pandas as pd
import json
import definitions


def extract_bioactivity_report(bioactivities_path):
    bioactivity_dict = dict()
    target_uniprot_dict = dict()
    for chunk in pd.read_csv(bioactivities_path, sep="\t", chunksize=10000):
        for index, row in chunk.iterrows():
            pdb = row["PDB"]
            row["CHAIN"]
            uniprot = row["UNIPROT"]
            row["INCHIKEY"]
            row["LIG"]
            mol_chembl_id = row["CHEMBL_ID_MOL"]
            tar_chembl_id = row["CHEMBL_ID_TAR"]
            if pdb not in bioactivity_dict:
                bioactivity_dict[pdb] = {tar_chembl_id: set([mol_chembl_id])}
                target_uniprot_dict[pdb] = {tar_chembl_id: set([uniprot])}
            else:
                if tar_chembl_id in bioactivity_dict[pdb]:
                    bioactivity_dict[pdb][tar_chembl_id].add(mol_chembl_id)
                    target_uniprot_dict[pdb][tar_chembl_id].add(uniprot)
                else:
                    bioactivity_dict[pdb][tar_chembl_id] = set([mol_chembl_id])
                    target_uniprot_dict[pdb][tar_chembl_id] = set([uniprot])
    return bioactivity_dict, target_uniprot_dict


def process_and_save_dicts(bioactivity_dict, target_uniprot_dict):
    for pdb in bioactivity_dict:
        for target in bioactivity_dict[pdb]:
            molecules = list(bioactivity_dict[pdb][target])
            uniprots = list(target_uniprot_dict[pdb][target])
            bioactivity_dict[pdb][target] = molecules
            target_uniprot_dict[pdb][target] = uniprots
    with open(definitions.output_files_path / "chembiograph.json", "w") as outfile:
        json.dump(bioactivity_dict, outfile)

    with open(
            definitions.output_files_path / "chembiograph_targets.json", "w"
    ) as outfile:
        json.dump(target_uniprot_dict, outfile)


if __name__ == "__main__":
    bioactivities_path = definitions.input_files_path / "bioactivity_report_full_28april.tsv"
    bioactivity_dict, target_uniprot_dict = extract_bioactivity_report(
        bioactivities_path
    )
    process_and_save_dicts(bioactivity_dict, target_uniprot_dict)
