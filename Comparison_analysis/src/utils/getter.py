import json

import pandas as pd

from src import definitions
from src.utils.db import query_or_readin_data, query_activitydb


def get_chembl_to_usmiles(recalc=False):
    chembl_mol_df = get_full_chembl_mol_df(recalc=recalc)
    chembl_to_usmiles = {}
    for index, row in chembl_mol_df.iterrows():
        chembl_to_usmiles[row["chembl_id"]] = row["canonical_smiles"]
    return chembl_to_usmiles

def get_activitydb_pdb_ligand_table(recalc=False):
    df_path = definitions.output_files_path / "activitydb_pdb_ligand_table.csv"
    sql_query = "SELECT * FROM pdb_ligand"
    chembl = False
    return query_or_readin_data(
        df_path=df_path, sql_query=sql_query, recalc=recalc, chembl=chembl
    )


def get_activitydb_file_processing(recalc=False):
    df_path = definitions.output_files_path / "pdbs_processed_activityfinder.csv"
    sql_query = "SELECT * FROM pdb_file_processing"
    chembl = False
    return query_or_readin_data(
        df_path=df_path, sql_query=sql_query, recalc=recalc, chembl=chembl
    )


def get_processed_pdbs(recalc=False):
    file_processing_df = get_activitydb_file_processing(recalc=recalc)
    processed_pdbs = set(file_processing_df["pdb_id"])
    return processed_pdbs


def get_valid_pdbs(recalc=False):
    file_processing_df = get_activitydb_file_processing(recalc=recalc)
    valid_pdbs = set(
        file_processing_df[file_processing_df["comment"] == "success"]["pdb_id"]
    )
    return valid_pdbs


def get_chembl_target_df(recalc=False):
    df_path = definitions.output_files_path / "chembl_target_dataframe.csv"
    sql_query = "SELECT * FROM target_dictionary"
    chembl = True
    return query_or_readin_data(
        df_path=df_path, sql_query=sql_query, recalc=recalc, chembl=chembl
    )


def get_valid_chembl_targets(recalc=False):
    target_df = get_chembl_target_df(recalc)
    valid_chembl_targets = [x for x in target_df["chembl_id"]]
    return valid_chembl_targets


def get_simple_chembl_mol_df(recalc=False):
    df_path = definitions.output_files_path / "valid_ligand_data.csv"
    sql_query = "SELECT * FROM molecule_dictionary"
    chembl = True
    return query_or_readin_data(
        df_path=df_path, sql_query=sql_query, recalc=recalc, chembl=chembl
    )


def get_full_chembl_mol_df(recalc=False):
    df_path = definitions.output_files_path / "molecule_dictionary.csv"
    sql_query = "SELECT molecule_dictionary.chembl_id, compound_structures.standard_inchi, compound_structures.canonical_smiles  FROM molecule_dictionary INNER JOIN compound_structures on molecule_dictionary.molregno = compound_structures.molregno"
    chembl = True
    return query_or_readin_data(
        df_path=df_path, sql_query=sql_query, recalc=recalc, chembl=chembl
    )


def get_valid_chembl_ligands(recalc=False):
    ligand_dataframe = get_simple_chembl_mol_df(recalc)
    return [x for x in ligand_dataframe["chembl_id"]]


def get_chembl_id_to_hetcode_dict(recalc=False):
    json_path = definitions.output_files_path / "biochem_chembl_to_hetcode.json"
    chembl_id_to_hetcode = {}
    if json_path.exists() and not recalc:
        with open(json_path, "r") as f:
            chembl_id_to_hetcode = json.load(f)
    else:
        chembio_data = pd.read_csv(
            definitions.biochemgraph_full, sep="\t"
        )
        for index, row in chembio_data.iterrows():
            chembl_id_mol = row["CHEMBL_ID_MOL"]
            hetcode = row["LIG"]
            if chembl_id_mol in chembl_id_to_hetcode.keys():
                continue
            else:
                chembl_id_to_hetcode[chembl_id_mol] = hetcode
    return chembl_id_to_hetcode


def get_chembl_ligands_in_activitydb(recalc=True):
    chembl_ligand_path = definitions.output_files_path / "chembl_ligands_activitydb.csv"
    if chembl_ligand_path.exists() and not recalc:
        chembl_ligand_data = pd.read_csv(chembl_ligand_path)
    else:
        sql_query = "SELECT * FROM chembl_ligand"
        chembl_ligand_data = query_activitydb(sql_query)
        chembl_ligand_data.to_csv(chembl_ligand_path)
    all_chembl_ligands = [x for x in chembl_ligand_data["chembl_id"]]
    return all_chembl_ligands

def get_chembl_targets_in_activitydb(recalc=True):
    chembl_target_path = definitions.output_files_path / "chembl_targets_activitydb.csv"
    if chembl_target_path.exists() and not recalc:
        chembl_target_data = pd.read_csv(chembl_target_path)
    else:
        sql_query = "SELECT * FROM target"
        chembl_target_data = query_activitydb(sql_query)
        chembl_target_data.to_csv(chembl_target_path)
    all_chembl_targets = [x for x in chembl_target_data["chembl_id"]]
    return all_chembl_targets


def get_sequence_identity_data_biochem():
    sequence_data = pd.read_csv(
        definitions.output_files_path / "sequence_identity_data_biochem_full.csv"
    )
    return sequence_data