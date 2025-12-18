import json

import pandas as pd
import plotnine as p9

from src import definitions

good_mutation_data = pd.read_csv(
    definitions.automated_path / "difference_variant_mutation_data_nogaps.csv"
)
print(good_mutation_data.shape)


n_pdbs = len(set([x for x in good_mutation_data["PDB"]]))
n_assays = len(set([x for x in good_mutation_data["assay"]]))
n_ligands = len(set([x[0:3] for x in good_mutation_data["name"]]))
print(n_pdbs, n_assays, n_ligands)
# Is the variant sequence mutationsfree? N_variant vs N_compoenent vs N_mutations
plot = p9.ggplot(good_mutation_data)
plot += p9.geom_histogram(p9.aes(x="n_variant"), binwidth=1, fill="red")
plot.save(definitions.plot_path / "Variant_sequence_mutations_variant.png")

plot = p9.ggplot(good_mutation_data)
plot += p9.geom_histogram(p9.aes(x="n_mutations"), binwidth=1, fill="green")
plot.save(definitions.plot_path / "Variant_sequence_mutations_mutations.png")

plot = p9.ggplot(good_mutation_data)
plot += p9.geom_histogram(p9.aes(x="n_component"), binwidth=1, fill="green")
plot.save(definitions.plot_path / "Variant_sequence_mutations_component.png")


general_data_dict = {}
for index, row in good_mutation_data.iterrows():
    pdb = row["PDB"]
    assay = row["assay"]
    name = str(row["name"])
    if pdb in general_data_dict.keys():
        if assay in general_data_dict[pdb]:
            if name not in general_data_dict[pdb][assay]:
                general_data_dict[pdb][assay] = {name: []}
        else:
            general_data_dict[pdb] = {assay: {name: []}}
    else:
        general_data_dict[pdb] = {assay: {name: []}}
print(general_data_dict)

for chunk in pd.read_csv(
    definitions.input_files_path / "parsed_activity_data_raw.csv", chunksize=10000
):
    for index, row in chunk.iterrows():
        pdb = row["pdb.pdb_id"]
        assay = row["assay.chembl_id"]
        ligand_hetcode = row["pdb_ligand.ligand_id"]
        ligand_id = row["chembl_ligand.chembl_id"]
        if pdb in general_data_dict:
            if assay in general_data_dict[pdb]:
                for current_name in general_data_dict[pdb][assay]:
                    if current_name[0:3] == ligand_hetcode:
                        if ligand_id not in general_data_dict[pdb][assay][current_name]:
                            general_data_dict[pdb][assay][current_name].append(
                                ligand_id
                            )
                            print(
                                f"Found ligand {ligand_id} for pdb {pdb} and assay {assay}"
                            )
                        break

with open(definitions.automated_path / "variant_sequence_data_dict.json", "w") as f:
    json.dump(general_data_dict, f)
