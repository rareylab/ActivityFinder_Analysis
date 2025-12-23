"""SQL queries for HET code variant analysis."""

# Query to compute distribution of unique InChI-Key variants per HET code
HET_CODE_DATA_EXTRACTION_SQL = """
/*
  data_extraction.sql

  Purpose:
    Compute the distribution of the number of variants per HET code (ligand_id),
    where a "variant" is now defined as a distinct raw NAOMI InChI-Key
    (column: standard_inchi_key) observed for that ligand_id in the pdb_ligand table.

  Output Columns:
    num_variants  - Number of distinct standard_inchi_key values for a ligand_id
    num_ligands   - How many ligand_ids exhibit exactly num_variants variants
    percentage    - (num_ligands / total number of ligand_ids) * 100

  Notes:
    - NULL standard_inchi_key values are ignored by COUNT(DISTINCT ...) (PostgreSQL behavior).
    - If desired, you can enforce exclusion explicitly by adding:
        WHERE standard_inchi_key IS NOT NULL
      inside the CTE.
    - Downstream plot/table code should annotate that the variant definition is
      "based on distinct raw NAOMI InChI-Keys (standard_inchi_key) per HET code".

  Revision History:
    - Updated to use standard_inchi_key instead of naomi_unique_smiles_chiral.
*/

WITH ligand_counts AS (
    SELECT
        ligand_id,
        COUNT(DISTINCT standard_inchi_key) AS unique_inchi_keys,
        COUNT(*) AS total_entries
    FROM
        pdb_ligand
    GROUP BY
        ligand_id
)
SELECT
    unique_inchi_keys AS num_variants,
    COUNT(*) AS num_ligands,
    CAST(COUNT(*) AS FLOAT)
        / (SELECT COUNT(*) FROM ligand_counts) * 100 AS percentage
FROM
    ligand_counts
GROUP BY
    unique_inchi_keys
ORDER BY
    unique_inchi_keys;
"""
