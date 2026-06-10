"""
fasta_utils.py
--------------
Utility functions for FASTA file processing and protein localization analysis.
Comparing DeepLoc outputs with Genome-scale Metabolic Model reaction databases

Dependencies:
    pandas, numpy, biopython

Usage:
    from fasta_utils import split_fasta, filter_hypothetical_proteins, count_sequences, rename_sequences
"""

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez


# ---------------------------------------------------------------------------
# FASTA Utilities
# ---------------------------------------------------------------------------

def batch_iterator(iterator, batch_size):
    """
    Yield successive batches from an iterator.

    Parameters
    ----------
    iterator : iterable
        Any iterable (e.g. SeqRecord objects from SeqIO.parse).
    batch_size : int
        Number of items per batch.

    Yields
    ------
    list
        A list of up to batch_size items.
    """
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def split_fasta(file_in, file_out, batch_size):
    """
    Split a large FASTA file into smaller files of a given sequence count.
    Useful when file is too large for NCBI BLAST

    Parameters
    ----------
    file_in : str
        Path to the input FASTA file.
    file_out : str
        Base name for output files (e.g. "output" → output_split_1.fasta).
    batch_size : int
        Number of sequences per output file.

    Returns
    -------
    int
        Number of split files created.

    Example
    -------
    split_fasta("proteins.fasta", "proteins_batch", 500)
    # Writes: proteins_batch_split_1.fasta, proteins_batch_split_2.fasta, ...
    """
    record_iter = SeqIO.parse(file_in, "fasta")
    num_files = 0

    for i, batch in enumerate(batch_iterator(record_iter, batch_size), start=1):
        new_file = f"{file_out}_split_{i}.fasta"
        with open(new_file, "w") as out_handle:
            count = SeqIO.write(batch, out_handle, "fasta")
        print(f"Wrote {count} records to {new_file}")
        num_files += 1

    return num_files


def filter_hypothetical_proteins(file_in, file_out, keyword="hypothetical protein"):
    """
    Extract sequences whose description contains a given keyword into a new FASTA file.

    Parameters
    ----------
    file_in : str
        Path to the input FASTA file.
    file_out : str
        Path for the output FASTA file containing matched sequences.
    keyword : str, optional
        Substring to search for in sequence descriptions (default: "hypothetical protein").

    Returns
    -------
    int
        Number of matching sequences written.

    Example
    -------
    filter_hypothetical_proteins("proteome.fasta", "hypothetical.fasta")
    """
    count = 0

    with open(file_in, "r") as infile, open(file_out, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if keyword.lower() in record.description.lower():
                outfile.write(f">{record.description}\n{record.seq}\n")
                count += 1

    print(f"{count} sequences matching '{keyword}' written to {file_out}")
    return count


def count_sequences(file_in):
    """
    Count the number of sequences in a FASTA file.

    Parameters
    ----------
    file_in : str
        Path to the input FASTA file.

    Returns
    -------
    int
        Number of sequences found.
    """
    count = 0
    with open(file_in, "r") as file:
        for line in file:
            if line.startswith(">"):
                count += 1

    print(f"{count} sequences in {file_in}")
    return count


def rename_sequences(file_in, file_out, prefix, dict_out="fasta_dictionary.txt"):
    """
    Rename FASTA sequence IDs to a standardised format: PREFIX + zero-padded index.
    Used for assigning GPRs to Metabolic model

    Parameters
    ----------
    file_in : str
        Path to the input FASTA file.
    file_out : str
        Path for the renamed output FASTA file.
    prefix : str
        Exactly 4-character prefix for new sequence IDs for specific organism (e.g. "PROT").
    dict_out : str, optional
        Path for the mapping file (original ID → new ID). Default: "fasta_dictionary.txt".

    Returns
    -------
    int
        Number of sequences renamed.

    Raises
    ------
    ValueError
        If prefix is not exactly 4 characters.

    Example
    -------
    rename_sequences("proteome.fasta", "proteome_renamed.fasta", "PROT")
    # >long_original_name  →  >PROT00001
    """
    if len(prefix) != 4:
        raise ValueError(f"Prefix must be exactly 4 characters, got '{prefix}' ({len(prefix)} chars).")

    count = 0
    with open(file_in, "r") as infile, \
         open(file_out, "w") as outfile, \
         open(dict_out, "w") as dictfile:

        for i, record in enumerate(SeqIO.parse(infile, "fasta"), start=1):
            new_id = f"{prefix}{i:05d}"
            outfile.write(f">{new_id}\n{record.seq}\n")
            dictfile.write(f"{record.description}\t{new_id}\n")
            count += 1

    print(f"Renamed {count} sequences. Dictionary saved to {dict_out}")
    return count


# ---------------------------------------------------------------------------
# Protein Localization Analysis
# ---------------------------------------------------------------------------

def load_localization_data(model_file, localization_file):
    """
    Load reaction and protein localization data from Excel files.

    Parameters
    ----------
    model_file : str
        Path to the metabolic model Excel file.
        Expected sheet: 'Reaction List' with columns: Abbreviation, Gene, GPR, Reaction Location.
    localization_file : str
        Path to the localization predictions Excel file.
        Expected sheet: 'Predictions' with columns: Protein ID, Predicted Localization.

    Returns
    -------
    tuple of (pd.DataFrame, pd.DataFrame)
        reactions_df, localization_df
    """
    reactions_df = pd.read_excel(
        model_file,
        sheet_name="Reaction List",
        usecols=["Abbreviation", "Gene", "GPR", "Reaction Location"]
    ).replace(np.nan, "", regex=True)

    localization_df = pd.read_excel(
        localization_file,
        sheet_name="Predictions",
        usecols=["Protein ID", "Predicted Localization"]
    )

    return reactions_df, localization_df


def map_genes_to_model(reactions_df, localization_df):
    """
    Flag which proteins from the localization file appear in the metabolic model,
    and count how many reactions each protein is associated with.

    Parameters
    ----------
    reactions_df : pd.DataFrame
        Reaction data from the metabolic model.
    localization_df : pd.DataFrame
        Protein localization predictions.

    Returns
    -------
    pd.DataFrame
        localization_df with two new columns:
        - 'In_Model': 'Y' if the protein appears in any reaction, else 'N'.
        - 'Reaction_Count': number of reactions the protein is associated with.
    """
    gene_in_model = []
    reaction_counts = []

    for protein_id in localization_df["Protein ID"]:
        count = reactions_df["Gene"].str.contains(protein_id, regex=False).sum()
        gene_in_model.append("Y" if count > 0 else "N")
        reaction_counts.append(count)

    localization_df = localization_df.copy()
    localization_df["In_Model"] = gene_in_model
    localization_df["Reaction_Count"] = reaction_counts

    print(f"Total Protein-Reaction matches: {sum(reaction_counts)}")
    return localization_df


def compare_localizations(reactions_df, localization_df):
    """
    Build a comparison table of predicted protein localization vs. reaction localization
    in the metabolic model, and flag mismatches.

    Parameters
    ----------
    reactions_df : pd.DataFrame
        Reaction data including 'Reaction Location' and 'Gene' columns.
    localization_df : pd.DataFrame
        Protein data including 'Protein ID' and 'Predicted Localization' columns.

    Returns
    -------
    tuple of (pd.DataFrame, pd.DataFrame)
        full_df    : all protein-reaction pairs with a 'Localization_Match' boolean column.
        mismatch_df: subset of full_df where localization does not match.
    """
    rows = []

    for _, loc_row in localization_df.iterrows():
        protein_id = loc_row["Protein ID"]
        protein_loc = loc_row["Predicted Localization"]

        matched_reactions = reactions_df[reactions_df["Gene"].str.contains(protein_id, regex=False)]

        for _, rxn_row in matched_reactions.iterrows():
            rows.append({
                "Protein_ID": protein_id,
                "Rxn_Abbreviation": rxn_row["Abbreviation"],
                "Reaction_Location": rxn_row["Reaction Location"],
                "Protein_Location": protein_loc,
                "Localization_Match": rxn_row["Reaction Location"] in protein_loc
            })

    full_df = pd.DataFrame(rows)
    mismatch_df = full_df[~full_df["Localization_Match"]].reset_index(drop=True)

    print(f"Total pairs: {len(full_df)} | Mismatches: {len(mismatch_df)}")
    return full_df, mismatch_df


# ---------------------------------------------------------------------------
# Example usage 
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    # --- FASTA processing ---
    split_fasta("proteome.fasta", "proteome_batch", batch_size=500)
    filter_hypothetical_proteins("proteome.fasta", "hypothetical.fasta")
    count_sequences("proteome.fasta")
    rename_sequences("proteome.fasta", "proteome_renamed.fasta", prefix="PROT")

    # --- Localization analysis ---
    reactions, localizations = load_localization_data(
        model_file="model.xlsx",
        localization_file="localization_predictions.xlsx"
    )
    localizations = map_genes_to_model(reactions, localizations)
    localizations.to_excel("gene_model_mapping.xlsx", index=False)

    full_comparison, mismatches = compare_localizations(reactions, localizations)
    mismatches.to_excel("localization_mismatches.xlsx", index=False)
