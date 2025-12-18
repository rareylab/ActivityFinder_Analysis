import definitions
import pandas as pd
import json

from src.utils import db
from mutation_analysis.extract_mutation_data import get_mutations_for_pdb
from mutation_analysis.analyse_mutation_data import read_mutation_data
from pathlib import Path


recalc = False
assays_with_variant_path = (
        definitions.output_files_path / "activitydb_assays_with_variant_id.csv"
)
if not assays_with_variant_path.exists() or recalc:
    # sql_query = "SELECT chembl_id, variant_id FROM assay where variant_id IS NOT NULL"
    sql_query = "SELECT a.chembl_id, a.variant_id, cvs.mutation FROM assay a JOIN chembl_variant_sequence cvs ON a.variant_id = cvs.cvsid WHERE a.variant_id IS NOT NULL"
    data = db.query_activitydb(sql_query)
    data.to_csv(assays_with_variant_path)
else:
    data = pd.read_csv(assays_with_variant_path)

relevant_assays = [x for x in data["chembl_id"]]
assay_to_pdb_data_path = definitions.output_files_path / "assay_to_data_dict.json"

if assay_to_pdb_data_path.exists() and not recalc:
    with open(assay_to_pdb_data_path, "r") as f:
        assay_ta_data_dict = json.load(f)
else:
    assay_ta_data_dict = {}
    for chunk in pd.read_csv(
            definitions.activities_raw,
        chunksize=10000,
        dtype=str,
    ):
        print("Chunki chunki chunk")
        for index, row in chunk.iterrows():
            assay = row["assay.chembl_id"]
            if not assay in relevant_assays:
                continue
            chembl_ligand = row["chembl_ligand.chembl_id"]
            pdb = row["pdb.pdb_id"]
            hetcode = row["pdb_ligand.ligand_id"]
            target = row["target.chembl_id"]
            if assay in assay_ta_data_dict.keys():
                if not pdb in assay_ta_data_dict[assay]:
                    assay_ta_data_dict[assay].append(pdb)
            else:
                assay_ta_data_dict[assay] = [pdb]
    with open(assay_to_pdb_data_path, "w") as f:
        json.dump(assay_ta_data_dict, f)

print(len(assay_ta_data_dict), len(relevant_assays))

variant_mutation_data_path = Path("/work/gutermuth/activityfinderanalysis/variant_data")
if not variant_mutation_data_path.exists():
    variant_mutation_data_path.mkdir()
different_mutation_data_path = (
        definitions.automated_path / "difference_variant_mutation_data.csv"
)

if not different_mutation_data_path.exists() or recalc:
    print("Start with mutation extraction")
    different_mutation_data = []
    for assay in assay_ta_data_dict.keys():
        print(assay)
        for pdb in assay_ta_data_dict[assay]:
            variant_path = variant_mutation_data_path / f"{pdb}_variant.csv"
            component_path = variant_mutation_data_path / f"{pdb}_component.csv"
            if not variant_path.exists() or not component_path.exists():
                component_data, variant_data = get_mutations_for_pdb(pdb)
                variant_data.to_csv(variant_path)
                component_data.to_csv(component_path)
            mutatian_dict_variant = read_mutation_data(variant_path)
            mutatian_dict_component = read_mutation_data(component_path)
            for ligand in mutatian_dict_component:
                for assay in mutatian_dict_component[ligand]:
                    if assay in relevant_assays:
                        for sequence in mutatian_dict_component[ligand][assay]:
                            n_mutations_component = len(
                                mutatian_dict_component[ligand][assay][sequence]
                            )
                            if (
                                ligand in mutatian_dict_variant
                                and assay in mutatian_dict_variant[ligand]
                                and sequence in mutatian_dict_variant[ligand][assay]
                            ):
                                n_mutations_variant = len(
                                    mutatian_dict_variant[ligand][assay][sequence]
                                )
                            else:
                                n_mutations_variant = 0
                            if n_mutations_variant < n_mutations_component:
                                print(
                                    pdb,
                                    ligand,
                                    assay,
                                    sequence,
                                    n_mutations_component,
                                    n_mutations_variant,
                                )
                                mutations = mutatian_dict_component[ligand][assay][
                                    sequence
                                ]
                                print(mutations)
                                mutations_variant = []
                                if not ligand in mutatian_dict_variant:
                                    print(f"Ligand {ligand} not in mutation data")
                                elif not assay in mutatian_dict_variant[ligand]:
                                    print(f"Assay {assay} not in mutation data")
                                elif (
                                    not sequence in mutatian_dict_variant[ligand][assay]
                                ):
                                    print(f"Sequence {sequence} not in mutation data")
                                else:
                                    mutations_variant = mutatian_dict_variant[ligand][
                                        assay
                                    ][sequence]
                                print(mutations_variant)
                                missing_mutations = [
                                    x
                                    for x in mutations
                                    if not x[0] in [y[0] for y in mutations_variant]
                                ]
                                if len(missing_mutations) > 0:
                                    different_mutation_data.append(
                                        [
                                            pdb,
                                            ligand,
                                            assay,
                                            sequence,
                                            n_mutations_component,
                                            n_mutations_variant,
                                            missing_mutations,
                                        ]
                                    )
                                    print(missing_mutations)
        different_mutation_data = pd.DataFrame(
            different_mutation_data,
            columns=[
                "PDB",
                "Ligand",
                "Assay",
                "Sequence",
                "N mutations component",
                "N mutations variant",
                "Missing mutations",
            ],
        )
        if different_mutation_data_path.exists():
            different_mutation_data.to_csv(
                different_mutation_data_path,
                mode="a",
                index=False,
                header=None
            )
        else:
            different_mutation_data.to_csv(different_mutation_data_path, index=False)
        different_mutation_data = []
else:
    print("Found existing mutation data")
    different_mutation_data = pd.read_csv(different_mutation_data_path)
