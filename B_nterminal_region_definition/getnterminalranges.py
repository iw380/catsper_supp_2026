import pandas as pd
import re
import os
import subprocess
from Bio import AlignIO, SeqIO

def get_n_terminal_slice_position(alignment_file):
    """
    Finds the start of the first highly conserved BLOCK in an alignment.
    This is a more robust way to define the end of the N-terminal region.
    """
    alignment = AlignIO.read(alignment_file, "fasta")
    num_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()
    conservation_threshold = 0.8
    min_consecutive_cols = 10
    consecutive_count = 0

    for i in range(alignment_length):
        column = alignment[:, i]
        non_gap_count = sum(1 for char in column if char != '-')
        if (non_gap_count / num_sequences) >= conservation_threshold:
            consecutive_count += 1
        else:
            consecutive_count = 0
            
        if consecutive_count >= min_consecutive_cols:
            core_start_position = i - min_consecutive_cols + 1
            print(f"Found start of conserved block at alignment position: {core_start_position}")
            return core_start_position
    
    print("Warning: Could not find a sufficiently long conserved block.")
    return 0

def calculate_nterminal_ranges(data_tsv, output_file):
    """
    Performs alignment and calculates the N-terminal residue range for each sequence.
    """
    # 1. Load the FULL-LENGTH sequences from your main data file
    print(f"Reading full-length sequences from '{data_tsv}'...")
    try:
        df = pd.read_csv(data_tsv, sep='\t')
    except FileNotFoundError:
        print(f"--- ERROR: Data file '{data_tsv}' not found. ---")
        return
        
    df_catsper1 = df[df['GeneName'] == 'CATSPER1'].copy().reset_index(drop=True)

    # 2. Create a temporary FASTA file for alignment
    print("\nCreating temporary FASTA file for alignment...")
    fasta_input_file = "temp_catsper1_full.fasta"
    with open(fasta_input_file, "w") as f:
        for index, row in df_catsper1.iterrows():
            clean_sequence = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', str(row['Sequence']).upper())
            if clean_sequence:
                f.write(f">{index}\n{clean_sequence}\n")

    # 3. Run MUSCLE alignment
    print("Running MUSCLE alignment...")
    aligned_fasta_file = "aligned_catsper1_full.fasta"
    try:
        subprocess.run(["muscle", "-align", fasta_input_file, "-output", aligned_fasta_file], check=True, capture_output=True, text=True)
        print("Alignment complete.")
    except Exception as e:
        print(f"\n--- MUSCLE ERROR ---: {e}")
        return

    # 4. Find the slice position from the alignment
    print("Analyzing alignment to define N-terminal regions...")
    slice_pos = get_n_terminal_slice_position(aligned_fasta_file)
    if slice_pos == 0:
        print("--- ANALYSIS FAILED: Could not determine a valid core start position. ---")
        return

    # 5. Calculate the N-terminal length for each species
    print("Calculating N-terminal length for each species...")
    alignment = {record.id: str(record.seq) for record in SeqIO.parse(aligned_fasta_file, "fasta")}
    
    results = []
    for index, row in df_catsper1.iterrows():
        if str(index) not in alignment:
            continue

        aligned_seq = alignment[str(index)]
        # Count non-gap characters up to the slice position to find the true length
        n_terminal_length = len(aligned_seq[:slice_pos].replace('-', ''))
        
        # Only add if the length is greater than 0
        if n_terminal_length > 0:
            results.append({
                'SpeciesName': row['SpeciesName'],
                'N_Terminal_Length': n_terminal_length,
                'ChimeraX_Selection': f"1-{n_terminal_length}"
            })

    # 6. Create a DataFrame, save to file, and print to screen
    results_df = pd.DataFrame(results)
    
    print("\n--- N-Terminal Selection Ranges for ChimeraX ---")
    print(results_df.to_string())

    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nResults also saved to '{output_file}'")

    # 7. Clean up temporary files
    os.remove(fasta_input_file)
    os.remove(aligned_fasta_file)
    print("Temporary files removed.")


if __name__ == '__main__':
    # --- Configuration ---
    # IMPORTANT: This script needs the original data file with the FULL sequences
    DATA_TSV_FILE = 'parsed_catsper_data.tsv' 
    OUTPUT_FILE = 'catsper1_nterminal_ranges.tsv'

    try:
        import pandas
        from Bio import AlignIO, SeqIO
    except ImportError as e:
        print(f"Error: Required library '{e.name}' is not installed.")
        print("Please install necessary libraries: pip install pandas biopython")
    else:
        calculate_nterminal_ranges(DATA_TSV_FILE, OUTPUT_FILE)