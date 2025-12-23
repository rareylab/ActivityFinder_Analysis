"""SQL queries for variant mutation analysis."""

# Query to analyze distribution of mutation counts in ChEMBL variant sequences
VARIANT_MUTATION_DATA_EXTRACTION_SQL = """
-- ChEMBL Variant Mutation Count Distribution Analysis
-- This query analyzes the distribution of mutation counts in ChEMBL variant sequences

WITH mutation_counts AS (
    SELECT
        accession,
        mutation,
        CASE
            WHEN mutation IS NULL OR mutation = '' THEN 0
            ELSE array_length(regexp_split_to_array(mutation, '[;,]'), 1)
        END AS mutation_count
    FROM chembl_33.variant_sequences
)
SELECT
    mutation_count AS num_mutations,
    COUNT(*) AS num_sequences,
    CAST(COUNT(*) AS FLOAT) / (
        SELECT COUNT(*) FROM mutation_counts
    ) * 100 AS percentage
FROM
    mutation_counts
GROUP BY
    mutation_count
ORDER BY
    mutation_count;
"""
