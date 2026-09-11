import pandas as pd
import re
import os
from Bio.SeqIO import parse
from Bio.PDB import PDBParser, SASA
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# --- ANALYSIS CONFIGURATION ---
RSA_THRESHOLD = 25.0 
MAX_ASA = {
    'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0, 'CYS': 167.0, 
    'GLN': 225.0, 'GLU': 223.0, 'GLY': 104.0, 'HIS': 224.0, 'HIE': 224.0, 
    'HID': 224.0, 'HIP': 224.0, 'ILE': 197.0, 'LEU': 201.0, 'LYS': 236.0, 
    'MET': 204.0, 'PHE': 240.0, 'PRO': 159.0, 'SER': 155.0, 'THR': 172.0, 
    'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0
}
# -----------------------------------------------------------

def process_full_length_data(data_tsv):
    """Gets total CATSPER1 length from the main parsed data file."""
    try:
        df = pd.read_csv(data_tsv, sep='\t')
        df_catsper1 = df[df['GeneName'] == 'CATSPER1'].copy()
        # In case of multiple entries, take the first/longest one
        df_agg = df_catsper1.loc[df_catsper1.groupby('SpeciesName')['SequenceLength'].idxmax()]
        return df_agg[['SpeciesName', 'SequenceLength']].rename(columns={'SequenceLength': 'Total_CATSPER1_Length'})
    except FileNotFoundError: return None

def process_nterminal_fasta(fasta_file):
    """Processes the N-terminal FASTA to get length, pI, and histidine percentage."""
    data_list = []
    try:
        for record in parse(fasta_file, "fasta"):
            sequence = str(record.seq).upper()
            if not sequence: continue
            species_name = record.id.rsplit('_', 1)[0].replace('_', ' ')
            n_term_len = len(sequence)
            his_percent = (sequence.count('H') / n_term_len) * 100 if n_term_len > 0 else 0
            
            try:
                pI = ProteinAnalysis(sequence.replace('X','')).isoelectric_point()
            except:
                pI = None # Handle sequences that might cause errors
                
            data_list.append({
                'SpeciesName': species_name, 'N_Terminal_Length': n_term_len,
                'N_Terminal_pI': pI, 'N_Terminal_Histidine_Percent': his_percent
            })
    except FileNotFoundError: return None
    return pd.DataFrame(data_list)

def process_surface_histidine_data(pdb_folder, manifest_file):
    """Processes all PDBs to get surface histidine counts."""
    analysis_results = []
    try:
        with open(manifest_file, 'r') as f:
            pdb_manifest = {line.strip().split(',', 1)[0].strip(): line.strip().split(',', 1)[1].strip() for line in f if ',' in line}
    except FileNotFoundError: return None

    for filename, species_name in pdb_manifest.items():
        pdb_path = os.path.join(pdb_folder, filename)
        if not os.path.exists(pdb_path): continue
        
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("protein", pdb_path)
            sr = SASA.ShrakeRupley()
            sr.compute(structure, level="R")

            surface_his_count, total_his_count = 0, 0
            for residue in structure[0].get_residues():
                if residue.get_resname() in {'HIS', 'HIE', 'HID', 'HIP'}:
                    total_his_count += 1
                    max_sasa = MAX_ASA.get(residue.get_resname(), 0)
                    if max_sasa > 0 and residue.sasa is not None and (residue.sasa / max_sasa * 100) > RSA_THRESHOLD:
                        surface_his_count += 1
            
            analysis_results.append({
                'SpeciesName': species_name,
                'Surface_Histidines': surface_his_count,
                'Total_Histidines': total_his_count
            })
        except Exception: continue
    return pd.DataFrame(analysis_results)

def parse_temperature_data(temp_file):
    """Parses the temperature data file for body temperature."""
    temp_list = []
    try:
        with open(temp_file, 'r') as f:
            next(f, None)
            for line in f:
                if ':' not in line: continue
                species_part, all_temps_part = line.split(':', 1)
                species_name = species_part.strip()
                body_temp_str = all_temps_part.split(',')[0].strip()
                def _parse_range(range_str):
                    if 'unknown' in range_str.lower(): return None
                    numbers = [float(n) for n in re.findall(r"[-+]?\d*\.\d+|\d+", range_str)]
                    return sum(numbers) / len(numbers) if numbers else None
                temp_list.append({'SpeciesName': species_name, 'AvgBodyTemp': _parse_range(body_temp_str)})
    except FileNotFoundError: return None
    return pd.DataFrame(temp_list)

def compile_all_data():
    """Main function to load, process, merge, and save all data."""
    print("--- Starting Master Data Compilation ---")
    
    # --- Define file paths ---
    FULL_DATA_TSV = 'parsed_catsper_data.tsv'
    NTERMINAL_FASTA = 'catsper1_n_terminals.fasta'
    PDB_FOLDER = 'n_terminal_pdbs'
    MANIFEST_FILE = 'pdb_manifest.txt'
    TEMP_FILE = 'species_temperatures.txt'
    
    # --- Load and process all data sources ---
    print("\n[1/4] Processing N-Terminal sequences...")
    n_terminal_df = process_nterminal_fasta(NTERMINAL_FASTA)
    
    print("\n[2/4] Processing full-length sequences...")
    full_length_df = process_full_length_data(FULL_DATA_TSV)
    
    print("\n[3/4] Processing temperature data...")
    temp_df = parse_temperature_data(TEMP_FILE)
    
    print("\n[4/4] Processing PDB files for surface histidine data...")
    surface_his_df = process_surface_histidine_data(PDB_FOLDER, MANIFEST_FILE)

    if n_terminal_df is None or full_length_df is None or temp_df is None:
        print("\n--- ERROR: One or more essential input files are missing. Aborting. ---")
        return

    # --- Merge all datasets ---
    print("\nMerging all data sources...")
    # Start with the N-terminal data as the base
    base_df = n_terminal_df
    
    # Function for case-insensitive merge
    def merge_data(df1, df2, on_col='SpeciesName'):
        df1['join_key'] = df1[on_col].str.lower()
        df2['join_key'] = df2[on_col].str.lower()
        merged = pd.merge(df1, df2.drop(columns=on_col), on='join_key', how='left')
        return merged.drop(columns='join_key')

    master_df = merge_data(base_df, full_length_df)
    master_df = merge_data(master_df, temp_df)
    if surface_his_df is not None:
        master_df = merge_data(master_df, surface_his_df)

    # --- Create and save the two final CSV files ---
    # CSV 1: The Primary Dataset
    primary_cols = [
        'SpeciesName', 'Total_CATSPER1_Length', 'N_Terminal_Length',
        'N_Terminal_pI', 'AvgBodyTemp', 'N_Terminal_Histidine_Percent'
    ]
    primary_df = master_df[primary_cols].dropna(subset=['N_Terminal_Length']).sort_values('SpeciesName').reset_index(drop=True)
    primary_output_file = 'master_dataset_primary.csv'
    primary_df.to_csv(primary_output_file, index=False)
    print(f"\n--- Primary Master Dataset --- (Saved to {primary_output_file})")
    print(primary_df.to_string())

    # CSV 2: The Surface Histidine Dataset
    if surface_his_df is not None:
        surface_cols = ['SpeciesName', 'Total_Histidines', 'Surface_Histidines', 'AvgBodyTemp']
        surface_df = master_df[surface_cols].dropna(subset=['Total_Histidines']).sort_values('SpeciesName').reset_index(drop=True)
        surface_output_file = 'master_dataset_surface_histidines.csv'
        surface_df.to_csv(surface_output_file, index=False)
        print(f"\n--- Surface Histidine Dataset --- (Saved to {surface_output_file})")
        print(surface_df.to_string())
    else:
        print("\nSkipping surface histidine dataset because PDB analysis data was not available.")

    print("\n--- Compilation Complete! ---")

if __name__ == '__main__':
    try:
        import pandas
        from Bio.SeqIO import parse
        from Bio.PDB import PDBParser, SASA
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except ImportError as e:
        print(f"Error: Required library '{e.name}' is not installed.")
        print("Please run: pip install pandas biopython")
    else:
        compile_all_data()