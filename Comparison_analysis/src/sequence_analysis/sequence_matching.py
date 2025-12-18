from Bio import Align
import json
import pandas as pd
import itertools
import requests

from src import definitions


def get_uniprot_for_chain(pdb_id, chain_id):
    new_results = requests.get(
        f"https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{pdb_id.lower()}/{chain_id}"
    )
    new_entity_id = new_results.json()[
        "rcsb_polymer_entity_instance_container_identifiers"
    ]["entity_id"]
    uniprot_id_raw = requests.get(
        f"https://data.rcsb.org/rest/v1/core/uniprot/{pdb_id.lower()}/{new_entity_id}"
    )
    uniprot_id = uniprot_id_raw.json()[0]["rcsb_uniprot_container_identifiers"][
        "uniprot_id"
    ]
    return uniprot_id


def get_identifier_from_uniprot(uniprot_id):
    result = requests.get(f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}")
    pfam_id = None
    for bla in result.json():
        if "Pfam" in result.json()[bla]:
            print(bla)
    for reference in result.json()["dbReferences"]:
        if reference["type"] == "Pfam":
            pfam_id = reference["id"]
    return pfam_id


def get_sequence_from_uniprot(uniprot_id):
    result = requests.get(f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}")
    data_json = result.json()
    return data_json["sequence"]["sequence"]


def construct_match(sequence1, sequence2):

    # Create a PairwiseAligner object
    aligner = Align.PairwiseAligner()

    # Set the alignment mode to global
    aligner.mode = "global"

    # Set the match and mismatch scores
    aligner.match_score = 2
    aligner.mismatch_score = -1

    # Set the affine gap costs
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    # Perform the alignment
    alignments = aligner.align(sequence1, sequence2)

    # Use the first alignment
    alignment = alignments[0]
    # print(alignment)

    # Extract aligned sequences
    aligned_query = alignment.aligned[0]
    aligned_subject = alignment.aligned[1]

    # Convert to string-based alignments
    aligned_query_seq = []
    aligned_subject_seq = []
    for (q_start, q_end), (s_start, s_end) in zip(aligned_query, aligned_subject):
        aligned_query_seq.append(sequence1[q_start:q_end])
        aligned_subject_seq.append(sequence2[s_start:s_end])

    aligned_query_str = "".join(aligned_query_seq)
    aligned_subject_str = "".join(aligned_subject_seq)

    sequence_identity = calculate_sequence_identity(
        aligned_query_str, aligned_subject_str
    )

    # Print the alignment and sequence identity
    # print("Aligned Query:   ", aligned_query_str)
    # print("Aligned Subject: ", aligned_subject_str)
    # print(f"Sequence Identity: {sequence_identity:.2f}%")
    return sequence_identity


# Calculate sequence identity for the query-covered region
def calculate_sequence_identity(aligned_query, aligned_subject):
    """
    Calculate the sequence identity between two aligned sequences.
    Sequence Identity is defined as:

        Sequence Identity (%) = (Number of Identical Residues / Number of Non-Gap Positions) * 100

    Where:
    - "Number of Identical Residues" is the count of matching residues between the aligned sequences.
    - "Number of Non-Gap Positions" is the total number of positions excluding gaps ('-') in the aligned sequences.

    This formula calculates the percentage of identical residues relative to the non-gap positions,
    making it independent of gaps introduced during alignment.

    Args:
        aligned_query:
            The aligned query sequence as a list of aligned residues.
        aligned_subject:
            The aligned hit sequence as a list of aligned residues.

    Returns:
        The sequence identity as a percentage value.
    """
    matches = 0
    length = 0
    for q, s in zip(aligned_query, aligned_subject):
        if q != "-" and s != "-":  # Ignore gaps
            length += 1
            if q == s:
                matches += 1
    if length == 0:  # Avoid division by zero
        return 0.0
    return matches / length * 100


def extract_seqres_data(seqres_data_path, relevant_pdbs):
    relevant_sequences = {}
    if definitions.extracted_seqres_data.exists():
        with open(definitions.extracted_seqres_data, "r") as f:
            relevant_sequences = json.load(f)
    else:
        with open(seqres_data_path, "r") as f:
            current_pdb = None
            current_chain = None
            current_sequence = None
            for line in f.readlines():
                if line.startswith(">"):
                    interesting_part = line.split()[0][1:]
                    pdb, chain = interesting_part.split("_")
                else:
                    current_sequence = str(line.strip())
                    if pdb in relevant_pdbs:
                        if pdb in relevant_sequences:
                            relevant_sequences[pdb][chain] = current_sequence
                        else:
                            relevant_sequences[pdb] = {chain: current_sequence}
        print(relevant_sequences)
        with open(definitions.extracted_seqres_data, "w") as f:
            json.dump(relevant_sequences, f)
    return relevant_sequences


# construct_match(pdb_sequence, uniprot_sequence)
def main():
    # Calculate the sequence identities for the complete chembio stuff
    biochem_input_path = definitions.biochemgraph_full
    seqres_path = definitions.input_files_path / "pdb_seqres_12aug.txt"
    biochem_data = pd.read_csv(biochem_input_path, sep="\t")
    print(biochem_data.head())
    combinations = set()
    uniprot_target_mapping = {}
    for index, row in biochem_data.iterrows():
        pdb = row["PDB"].upper()
        uniprot = row["UNIPROT"]
        current = pdb + "-" + uniprot
        target_id = row["CHEMBL_ID_TAR"]
        if uniprot in uniprot_target_mapping.keys():
            uniprot_target_mapping[uniprot].append(target_id)
        else:
            uniprot_target_mapping[uniprot] = [target_id]
        print(pdb, uniprot, current)
        combinations.add(current)
    print(len(combinations), combinations)
    relevant_pdbs = [x.split("-")[0].lower() for x in combinations]
    print(relevant_pdbs)
    relevant_sequences = extract_seqres_data(seqres_path, relevant_pdbs)

    sequence_identity_data = []
    for combination in combinations:
        pdb, uniprot = combination.split("-")
        if pdb.lower() not in relevant_sequences:
            print(f"Combination {combination} was skipped")
            continue
        uniprot_sequences = [get_sequence_from_uniprot(uniprot)]
        chains = [x for x in relevant_sequences[pdb.lower()]]
        pdb_sequences = [relevant_sequences[pdb.lower()][x] for x in chains]
        sequence_combinations = list(
            itertools.product(uniprot_sequences, pdb_sequences)
        )
        best_sequence_identity = 0
        best_pdb_sequence = None
        best_uniprot_sequence = None
        for sequence_comb in sequence_combinations:
            current_uniprot_sequence = sequence_comb[0].strip()
            current_pdb_sequence = sequence_comb[1].strip()
            current_sequence_identity = construct_match(
                current_uniprot_sequence, current_pdb_sequence
            )
            if current_sequence_identity > best_sequence_identity:
                best_sequence_identity = current_sequence_identity
                best_uniprot_sequence = sequence_comb[0]
                best_pdb_sequence = sequence_comb[1]
        for target_id in uniprot_target_mapping[uniprot]:
            sequence_identity_data.append(
                [
                    pdb,
                    uniprot,
                    target_id,
                    best_sequence_identity,
                    best_pdb_sequence,
                    best_uniprot_sequence,
                ]
            )

    sequence_identity_data = pd.DataFrame(
        sequence_identity_data,
        columns=[
            "PDB",
            "Uniprot ID",
            "ChEMBL Target",
            "Best sequence identity",
            "Best PDB seqres sequence",
            "Best uniprot sequence",
        ],
    )
    sequence_identity_data.drop_duplicates(inplace=True)
    sequence_identity_data.to_csv(
        definitions.output_files_path / "sequence_identity_data_biochem_full.csv"
    )


if __name__ == "__main__":
    main()
