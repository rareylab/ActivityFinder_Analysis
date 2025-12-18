from src import definitions
import pandas as pd
import ast

new_data = []
maximum = 6856140
n_header = 0
n_skipped_mutations = 0
n_skipped_rows = 0
for chunk in pd.read_csv(
        definitions.automated_path / "difference_variant_mutation_data.csv", chunksize=10000
):
    for index, row in chunk.iterrows():
        current_percentage = round(100 * index / maximum, 2)
        pdb = row["PDB"]
        ligand = row["Ligand"]
        assay = row["Assay"]
        sequence = row["Sequence"]
        n_component = row["N mutations component"]
        n_variant = row["N mutations variant"]
        raw_mutations = row["Missing mutations"]
        if raw_mutations == "Missing mutations":
            n_header += 1
            continue
        missing_mutations = ast.literal_eval(raw_mutations)
        mutations_to_drop = []
        for index, mutation in enumerate(missing_mutations):
            if mutation[1] == "_" or mutation[2] == "_":
                mutations_to_drop.append(index)
                n_skipped_mutations += 1
        filtered_mutations = [
            missing_mutations[i]
            for i in range(len(missing_mutations))
            if i not in mutations_to_drop
        ]
        mutated_residues = [x[0] for x in filtered_mutations]
        print(
            f"Currently at {current_percentage}% with {len(missing_mutations)} raw and {len(filtered_mutations)} left mutations and a difference of {len(missing_mutations) - len(filtered_mutations)}"
        )
        if len(filtered_mutations) >= 1:
            new_data.append(
                [
                    pdb,
                    ligand,
                    assay,
                    sequence,
                    missing_mutations,
                    mutated_residues,
                    len(missing_mutations),
                    n_component,
                    n_variant,
                ]
            )
        else:
            n_skipped_rows += 1

print( 
    f"There were {n_header} rows with a header, {n_skipped_mutations} mutations skipped which resulted in {n_skipped_rows} rows being skipped completely"
)
print(f"There are {len(new_data)} rows left")
new_data = pd.DataFrame(
    new_data,
    columns=[
        "PDB",
        "name",
        "assay",
        "sequence",
        "mutations",
        "residues",
        "n_mutations",
        "n_component",
        "n_variant",
    ],
)
new_data.to_csv(
    definitions.automated_path / "difference_variant_mutation_data_nogaps.csv"
)
