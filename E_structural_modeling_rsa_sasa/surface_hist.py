import pandas as pd
import plotly.express as px
import re
import os
from Bio.PDB import PDBParser, SASA

# --- ANALYSIS CONFIGURATION ---
RSA_THRESHOLD = 25.0 

# --- Maximum SASA values (from Tien et al., 2013) ---
MAX_ASA = {
    'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0,
    'CYS': 167.0, 'GLN': 225.0, 'GLU': 223.0, 'GLY': 104.0,
    'HIS': 224.0, 'HIE': 224.0, 'HID': 224.0, 'HIP': 224.0,
    'ILE': 197.0, 'LEU': 201.0, 'LYS': 236.0, 'MET': 204.0,
    'PHE': 240.0, 'PRO': 159.0, 'SER': 155.0, 'THR': 172.0,
    'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0
}
# -----------------------------------------------------------

def parse_dual_temperature_data(temp_file):
    """Parses the dual-temp file and returns a pandas DataFrame."""
    temp_list = []
    print(f"Reading dual temperature data from '{temp_file}'...")
    try:
        with open(temp_file, 'r') as f:
            next(f, None)
            for line in f:
                if ':' not in line: continue
                species_part, all_temps_part = line.split(':', 1)
                species_name = species_part.strip()
                temp_ranges = all_temps_part.split(',')
                body_temp_str = temp_ranges[0].strip()
                def _parse_range(range_str):
                    if 'unknown' in range_str.lower(): return None
                    numbers = [float(n) for n in re.findall(r"[-+]?\d*\.\d+|\d+", range_str)]
                    return sum(numbers) / len(numbers) if numbers else None
                temp_list.append({'SpeciesName': species_name, 'AvgBodyTemp': _parse_range(body_temp_str)})
    except FileNotFoundError: return None
    print(f"Successfully parsed temperatures for {len(temp_list)} species.")
    return pd.DataFrame(temp_list)

def analyze_pdb_file_with_rsa(pdb_path):
    """
    Uses Biopython's SASA calculator and a reference table to calculate RSA,
    then counts surface histidines.
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_path)
        sr = SASA.ShrakeRupley()
        sr.compute(structure, level="R")

        surface_histidine_count = 0
        total_histidine_count = 0 # <-- We will now also return this value
        total_residue_count = 0

        for residue in structure[0].get_residues():
            total_residue_count += 1
            res_name = residue.get_resname()
            
            if res_name in {'HIS', 'HIE', 'HID', 'HIP'}:
                total_histidine_count += 1 # <-- Count every histidine
                absolute_sasa = residue.sasa
                max_sasa = MAX_ASA.get(res_name, 0)
                
                if max_sasa > 0 and absolute_sasa is not None:
                    rsa = (absolute_sasa / max_sasa) * 100
                    if rsa > RSA_THRESHOLD:
                        surface_histidine_count += 1

        return {
            'Surface_Histidines': surface_histidine_count,
            'Total_Histidines': total_histidine_count, # <-- Add to output
            'Total_Residues': total_residue_count
        }

    except Exception as e:
        print(f"  - ERROR: Biopython failed to process {os.path.basename(pdb_path)}. Skipping. Error: {e}")
        return None

def read_pdb_manifest(manifest_file):
    """Reads the manifest file mapping PDB filenames to species names."""
    manifest = {}
    print(f"Reading PDB manifest from '{manifest_file}'...")
    try:
        with open(manifest_file, 'r') as f:
            for line in f:
                if ',' in line:
                    filename, species_name = line.strip().split(',', 1)
                    manifest[filename.strip()] = species_name.strip()
    except FileNotFoundError: return None
    print(f"Successfully loaded manifest for {len(manifest)} files.")
    return manifest

def create_surface_his_density_plot(pdb_folder, manifest_file, temp_file, output_html):
    """Main function to orchestrate the analysis and plotting."""
    pdb_manifest = read_pdb_manifest(manifest_file)
    temp_df = parse_dual_temperature_data(temp_file)
    if pdb_manifest is None or temp_df is None: return

    analysis_results = []
    print(f"\nAnalyzing PDB files based on manifest using Biopython (RSA method)...")
    
    for filename, species_name in pdb_manifest.items():
        pdb_path = os.path.join(pdb_folder, filename)
        if not os.path.exists(pdb_path): continue
        
        print(f"  -> Processing {filename} for species: {species_name}")
        analysis_data = analyze_pdb_file_with_rsa(pdb_path)
        
        if analysis_data:
            surface_his_density = (analysis_data['Surface_Histidines'] / analysis_data['Total_Residues'] * 100) if analysis_data['Total_Residues'] > 0 else 0
            analysis_results.append({
                'SpeciesName': species_name,
                'Total_Histidines': analysis_data['Total_Histidines'], # <-- Capture this value
                'Surface_Histidines': analysis_data['Surface_Histidines'],
                'Total_N_Terminal_Residues': analysis_data['Total_Residues'],
                'Surface_His_Density_RSA': surface_his_density
            })
    
    results_df = pd.DataFrame(analysis_results)
    results_df['join_key'] = results_df['SpeciesName'].str.lower()
    temp_df['join_key'] = temp_df['SpeciesName'].str.lower()
    final_df = pd.merge(results_df, temp_df.drop('SpeciesName', axis=1), on='join_key', how='inner').drop('join_key', axis=1)

    if final_df.empty: print("\n--- ERROR: Merged data is empty. ---"); return

    # --- UPDATED PRINTING SECTION ---
    # Select and reorder columns for a clearer console output
    print_columns = [
        'SpeciesName', 
        'AvgBodyTemp',
        'Total_Histidines', 
        'Surface_Histidines',
        'Surface_His_Density_RSA',
        'Total_N_Terminal_Residues'
    ]
    print(f"\nSuccessfully merged data for {len(final_df)} species:")
    print(final_df[print_columns].to_string())
    # ------------------------------

    df_plot = final_df.dropna(subset=['AvgBodyTemp'])
    fig = px.scatter(
        df_plot, x='AvgBodyTemp', y='Surface_His_Density_RSA', text='SpeciesName',
        hover_data=['SpeciesName', 'AvgBodyTemp', 'Surface_His_Density_RSA', 'Surface_Histidines', 'Total_Histidines', 'Total_N_Terminal_Residues'],
        title='<b>N-Terminal Surface Histidine Density (RSA > 25%) vs. Body Temperature</b>',
        labels={'AvgBodyTemp': 'Average Body Temperature (°C)', 'Surface_His_Density_RSA': 'Surface Histidine % (out of Total N-Term Residues)'}
    )
    fig.update_traces(textposition='top center', textfont_size=10)
    fig.update_layout(height=800, width=1200)
    
    fig.write_html(output_html)
    print(f"\nPlot saved to '{output_html}'")
    fig.show()

if __name__ == '__main__':
    PDB_FOLDER = 'n_terminal_pdbs'
    MANIFEST_FILE = 'pdb_manifest.txt'
    TEMP_FILE = 'species_temperatures.txt'
    OUTPUT_PLOT_HTML = 'surface_his_density_rsa_plot.html'

    try:
        import pandas, plotly
        from Bio.PDB import PDBParser, SASA
    except ImportError as e:
        print(f"Error: Required library '{e.name}' is not installed.")
        print("Please install necessary libraries: pip install pandas plotly biopython")
    else:
        create_surface_his_density_plot(PDB_FOLDER, MANIFEST_FILE, TEMP_FILE, OUTPUT_PLOT_HTML)