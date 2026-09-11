import pandas as pd
import re
import os
import subprocess
from Bio import AlignIO, SeqIO

# --- MANUAL OVERRIDE CONFIGURATION ---
# Set this to a number (e.g., 500) to manually define the alignment slice position.
# Set it to 'None' to use the automatic detection method.
MANUAL_SLICE_POSITION = None 
# ------------------------------------

def get_n_terminal_slice_position(alignment_file):
    """
    Finds the start of the first highly conserved BLOCK in an alignment.
    This is a more robust way to define the end of the N-terminal region.
    """
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
    except ValueError:
        print(f"--- ERROR: Alignment file '{alignment_file}' is empty or not in FASTA format. ---")
        return 0

    num_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()
    
    if num_sequences == 0 or alignment_length == 0:
        print("--- ERROR: Alignment is empty. ---")
        return 0

    # --- Robust "Core-Finding" Logic ---
    conservation_threshold = 0.8  # A column is "conserved" if >80% is not a gap.
    min_consecutive_cols = 10     # We need 10 conserved columns in a row to call it the "core".
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
    
    print("Warning: Could not find a sufficiently long conserved block. Using fallback.")
    return alignment_length // 2

def extract_n_terminals_to_fasta(data_tsv, output_fasta):
    """
    Performs alignment, slices the N-terminal regions, and saves them to a FASTA file.
    Allows for a manual override of the slice position.
    """
    # 1. Load sequence data and filter for CATSPER1
    try:
        df = pd.read_csv(data_tsv, sep='\t')
    except FileNotFoundError:
        print(f"--- ERROR: Data file '{data_tsv}' not found. ---")
        return
        
    df_catsper1 = df[df['GeneName'] == 'CATSPER1'].copy().reset_index()

    # 2. Create a temporary FASTA file for alignment with cleaned sequences
    print("\nCreating temporary FASTA file for alignment...")
    fasta_input_file = "temp_catsper1.fasta"
    with open(fasta_input_file, "w") as f:
        for index, row in df_catsper1.iterrows():
            clean_sequence = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', str(row['Sequence']).upper())
            if clean_sequence:
                f.write(f">{index}\n{clean_sequence}\n")

    # 3. Run MUSCLE alignment
    print("Running MUSCLE alignment...")
    aligned_fasta_file = "aligned_catsper1.fasta"
    try:
        # Using the modern MUSCLE v5 command-line options
        subprocess.run(["muscle", "-align", fasta_input_file, "-output", aligned_fasta_file], check=True, capture_output=True, text=True)
        print("Alignment complete.")
    except Exception as e:
        print(f"\n--- MUSCLE ERROR ---")
        print("An error occurred while running MUSCLE. Please ensure it is installed and in your PATH.")
        print(f"Error details: {e}")
        return

    # 4. Find the slice position from the alignment (or use the manual override)
    slice_pos = 0
    if MANUAL_SLICE_POSITION is not None and isinstance(MANUAL_SLICE_POSITION, int):
        slice_pos = MANUAL_SLICE_POSITION
        print(f"\nUsing MANUAL slice position: {slice_pos}")
    else:
        print("\nAnalyzing alignment to automatically define N-terminal regions...")
        slice_pos = get_n_terminal_slice_position(aligned_fasta_file)
        if slice_pos == 0:
            print("--- ANALYSIS FAILED: Could not determine a valid core start. Aborting. ---")
            return
    
    # 5. Extract N-terminals and save to the final output FASTA file
    print(f"Extracting N-terminal sequences and saving to '{output_fasta}'...")
    alignment = {record.id: str(record.seq) for record in SeqIO.parse(aligned_fasta_file, "fasta")}
    
    count_saved = 0
    with open(output_fasta, 'w') as f_out:
        for index, row in df_catsper1.iterrows():
            if str(index) not in alignment:
                continue

            aligned_seq = alignment[str(index)]
            # Count non-gap characters up to the slice position to find original length
            n_terminal_length = len(aligned_seq[:slice_pos].replace('-', ''))
            
            # Slice the *original* sequence
            n_terminal_seq = row['Sequence'][:n_terminal_length]
            
            if n_terminal_seq:
                species_name_formatted = row['SpeciesName'].replace(' ', '_')
                header = f">{species_name_formatted}_{row['index']}" # Use original index for consistency
                f_out.write(f"{header}\n{n_terminal_seq}\n")
                count_saved += 1

    print(f"Successfully saved {count_saved} N-terminal sequences.")

    # 6. Clean up temporary files
    os.remove(fasta_input_file)
    os.remove(aligned_fasta_file)
    print("Temporary files removed.")


if __name__ == '__main__':
    # --- Configuration ---
    DATA_TSV_FILE = 'parsed_catsper_data.tsv'
    OUTPUT_FASTA_FILE = 'catsper1_n_terminals.fasta'

    try:
        import pandas
        from Bio import AlignIO, SeqIO
    except ImportError as e:
        print(f"Error: Required library '{e.name}' is not installed.")
        print("Please install necessary libraries: pip install pandas biopython")
    else:
        extract_n_terminals_to_fasta(DATA_TSV_FILE, OUTPUT_FASTA_FILE)