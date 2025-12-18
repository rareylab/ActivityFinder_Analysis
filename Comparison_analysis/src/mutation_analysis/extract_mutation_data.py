import json
import progressbar
from src import definitions
from src.utils import db


def get_mutations_for_pdb(pdb_id):
    pdb_id = pdb_id.upper()
    sql_query_component = f"""SELECT DISTINCT blast_match.bmid, blast_match.query_sid, blast_match.hit_sid, 
    blast_match.percent_identity, blast_match.match_quality_score, pdb.pdb_id, active_site.naomi_ligand_id, 
    chain.name, sequence_type.type, target.chembl_id, assay.chembl_id, residue.residue_name, residue.insertion_code, 
    residue.sequence_number, quadruple.qid, quadruple.atom_seq_position, quadruple.atom_seq_residue, 
    quadruple.primary_seq_position, quadruple.primary_seq_residue, quadruple.blast_query_seq_position, 
    quadruple.blast_query_seq_residue, quadruple.blast_hit_seq_position, quadruple.blast_hit_seq_residue,
     'Component Sequence' AS chembl_target_sequence_source, assay.dsid AS assay_dsid, pdb.dsid AS pdb_dsid, 
     chembl_component_sequence.component_id AS chembl_sequence_id FROM blast_match 
     INNER JOIN aligned_sequence_blast_match ON blast_match.bmid = aligned_sequence_blast_match.bmid 
     INNER JOIN pdb_sequence_aligned_sequence_blast_match ON aligned_sequence_blast_match.asbmid = pdb_sequence_aligned_sequence_blast_match.asbmid 
     INNER JOIN pdb_sequence ON pdb_sequence_aligned_sequence_blast_match.psid = pdb_sequence.psid 
     INNER JOIN pdb ON pdb_sequence.pid = pdb.pid 
     INNER JOIN pdb_ligand_pdb_id ON pdb.pid = pdb_ligand_pdb_id.pid 
     INNER JOIN pdb_ligand ON pdb_ligand_pdb_id.pdblid = pdb_ligand.pdblid 
     INNER JOIN sequence_type ON pdb_sequence.stid = sequence_type.stid 
     INNER JOIN chain ON pdb_sequence.cid = chain.cid 
     INNER JOIN active_site ON active_site.pid = pdb.pid AND active_site.pdblid = pdb_ligand.pdblid 
     INNER JOIN active_site_residue_sequence ON active_site.actsid = active_site_residue_sequence.actsid 
     INNER JOIN residue ON residue.rid = active_site_residue_sequence.rid AND chain.cid = residue.cid 
     INNER JOIN quadruple ON aligned_sequence_blast_match.asbmid = quadruple.asbmid AND residue.residue_name = quadruple.atom_seq_residue AND residue.sequence_number = quadruple.atom_seq_res_seq_num 
     INNER JOIN chembl_component_sequence ON blast_match.hit_sid = chembl_component_sequence.sid 
     INNER JOIN sequence_target ON blast_match.hit_sid = sequence_target.sid 
     INNER JOIN target ON sequence_target.chembl_tid = target.chembl_tid 
     INNER JOIN assay ON target.tid = assay.tid 
     INNER JOIN activity ON assay.asid = activity.asid 
     INNER JOIN chembl_ligand ON activity.molregno = chembl_ligand.molregno 
     INNER JOIN pdb_ligand_chembl_ligand ON chembl_ligand.clid = pdb_ligand_chembl_ligand.clid 
     WHERE blast_match.percent_identity > 80 AND pdb.pdb_id = '{pdb_id}'"""
    mutation_data_component = db.query_activitydb(sql_query_component)
    sql_query_variant = (
        f"SELECT DISTINCT blast_match.bmid, blast_match.query_sid, blast_match.hit_sid, "
        f"blast_match.percent_identity, blast_match.match_quality_score, pdb.pdb_id, "
        f"active_site.naomi_ligand_id, chain.name, sequence_type.type, target.chembl_id, "
        f"assay.chembl_id, residue.residue_name, residue.insertion_code, residue.sequence_number, "
        f"quadruple.qid, quadruple.atom_seq_position, quadruple.atom_seq_residue, "
        f"quadruple.primary_seq_position, quadruple.primary_seq_residue, "
        f"quadruple.blast_query_seq_position, quadruple.blast_query_seq_residue, "
        f"quadruple.blast_hit_seq_position, quadruple.blast_hit_seq_residue, "
        f"'Variant Sequence' AS chembl_target_sequence_source, assay.dsid AS assay_dsid, "
        f"pdb.dsid AS pdb_dsid, chembl_variant_sequence.variant_id AS chembl_sequence_id FROM blast_match "
        f"INNER JOIN aligned_sequence_blast_match ON blast_match.bmid = aligned_sequence_blast_match.bmid "
        f"INNER JOIN pdb_sequence_aligned_sequence_blast_match ON aligned_sequence_blast_match.asbmid = pdb_sequence_aligned_sequence_blast_match.asbmid "
        f"INNER JOIN pdb_sequence ON pdb_sequence_aligned_sequence_blast_match.psid = pdb_sequence.psid INNER JOIN pdb ON pdb_sequence.pid = pdb.pid "
        f"INNER JOIN pdb_ligand_pdb_id ON pdb.pid = pdb_ligand_pdb_id.pid "
        f"INNER JOIN pdb_ligand ON pdb_ligand_pdb_id.pdblid = pdb_ligand.pdblid "
        f"INNER JOIN sequence_type ON pdb_sequence.stid = sequence_type.stid "
        f"INNER JOIN chain ON pdb_sequence.cid = chain.cid "
        f"INNER JOIN active_site ON active_site.pid = pdb.pid AND active_site.pdblid = pdb_ligand.pdblid "
        f"INNER JOIN active_site_residue_sequence ON active_site.actsid = active_site_residue_sequence.actsid "
        f"INNER JOIN residue ON residue.rid = active_site_residue_sequence.rid AND chain.cid = residue.cid "
        f"INNER JOIN quadruple ON aligned_sequence_blast_match.asbmid = quadruple.asbmid AND residue.residue_name = quadruple.atom_seq_residue AND residue.sequence_number = quadruple.atom_seq_res_seq_num "
        f"INNER JOIN chembl_variant_sequence ON blast_match.hit_sid = chembl_variant_sequence.sid "
        f"INNER JOIN sequence_target ON blast_match.hit_sid = sequence_target.sid "
        f"INNER JOIN target ON sequence_target.chembl_tid = target.chembl_tid "
        f"INNER JOIN assay ON chembl_variant_sequence.variant_id = assay.variant_id AND target.tid = assay.tid "
        f"INNER JOIN activity ON assay.asid = activity.asid "
        f"INNER JOIN chembl_ligand ON activity.molregno = chembl_ligand.molregno "
        f"INNER JOIN pdb_ligand_chembl_ligand ON chembl_ligand.clid = pdb_ligand_chembl_ligand.clid "
        f"WHERE blast_match.percent_identity > 80 AND pdb.pdb_id = '{pdb_id}'"
    )
    mutation_data_variant = db.query_activitydb(sql_query_variant)
    return mutation_data_component, mutation_data_variant


if "__main__" == __name__:
    with open(definitions.output_files_path / "chembiograph.json", "r") as f:
        chembio_full = json.load(f)

    if not definitions.mutation_data_path.exists():
        definitions.mutation_data_path.mkdir()

    all_pdbs = [x for x in chembio_full]
    n_pdbs = len(all_pdbs)
    for pdb in progressbar.progressbar(all_pdbs, max_value=n_pdbs):
        component_path = definitions.mutation_data_path / f"{pdb}_component.csv"
        variant_path = definitions.mutation_data_path / f"{pdb}_variant.csv"
        # if variant_path.exists() and component_path.exists():
        #    continue
        component, variant = get_mutations_for_pdb(pdb)
        component.to_csv(component_path)
        variant.to_csv(variant_path)
