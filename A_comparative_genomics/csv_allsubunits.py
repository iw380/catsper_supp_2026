import pandas as pd
import os

def create_cleaned_length_summary(input_tsv, species_to_include, output_csv):
    """
    Reads the main parsed data file, removes length outliers for each subunit,
    filters for a specific list of species, and writes the lengths to a
    wide-format CSV file.
    """
    # --- 1. Load the Data ---
    print(f"Reading sequence data from '{input_tsv}'...")
    try:
        all_seq_df = pd.read_csv(input_tsv, sep='\t')
    except FileNotFoundError:
        print(f"--- ERROR: Input file '{input_tsv}' not found. ---")
        return

    # --- 2. Remove Outliers using IQR Method (per subunit) ---
    print("\nRemoving length outliers for each subunit...")
    cleaned_dfs = []
    # Group the dataframe by the 'GeneName' column
    for gene, group_df in all_seq_df.groupby('GeneName'):
        Q1 = group_df['SequenceLength'].quantile(0.25)
        Q3 = group_df['SequenceLength'].quantile(0.75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR
        
        # Keep only the rows that are within the statistical fences
        inliers = group_df[(group_df['SequenceLength'] >= lower_bound) & (group_df['SequenceLength'] <= upper_bound)]
        outliers_removed = len(group_df) - len(inliers)
        if outliers_removed > 0:
            print(f"  - {gene}: Removed {outliers_removed} outlier(s).")
        cleaned_dfs.append(inliers)

    # Combine the cleaned dataframes back into one
    df_cleaned = pd.concat(cleaned_dfs)
    print(f"\nTotal entries before cleaning: {len(all_seq_df)}. After cleaning: {len(df_cleaned)}.")

    # --- 3. Filter by the Provided Species List ---
    print(f"\nFiltering for {len(species_to_include)} specified species...")
    df_filtered = df_cleaned[df_cleaned['SpeciesName'].isin(species_to_include)].copy()

    if df_filtered.empty:
        print("\n--- WARNING: No matching species were found in the cleaned data. ---")
        print("Please check the spelling in your SPECIES_TO_INCLUDE list.")
        return
        
    print(f"Found data for {len(df_filtered['SpeciesName'].unique())} of the specified species.")

    # --- 4. Reshape the Data to Wide Format ---
    df_pivot = df_filtered.pivot_table(
        index='SpeciesName', 
        columns='GeneName', 
        values='SequenceLength',
        aggfunc='mean' 
    )
    df_pivot.reset_index(inplace=True)
    df_pivot.columns.name = None

    # --- 5. Sort columns for readability ---
    other_cols = sorted([col for col in df_pivot.columns if col != 'SpeciesName'])
    sorted_cols = ['SpeciesName'] + other_cols
    df_pivot = df_pivot[sorted_cols]

    # --- 6. Save the Final CSV File ---
    try:
        df_pivot.to_csv(output_csv, index=False)
        print(f"\nSuccessfully saved the cleaned summary to '{output_csv}'.")
        print("\n--- Data Preview ---")
        print(df_pivot.to_string())
    except Exception as e:
        print(f"\n--- ERROR: Could not save the file. Error: {e}")


if __name__ == '__main__':
    # --- CONFIGURATION ---
    # 1. Define the list of species you want to include in your final output.
    SPECIES_TO_INCLUDE = [
        "Acropora digitifera",
    "Acropora millepora",
    "Alligator mississippiensis",
    "Apostichopus japonicus",
    "Bos taurus",
    "Canis lupus familiaris",
    "Capra hircus",
    "Capricornis sumatraensis",
    "Cavia porcellus",
    "Ciona intestinalis",
    "Elephantulus edwardii",
    "Equus caballus",
    "Fukomys damarensis",
    "Heterocephalus glaber",
    "Homo sapiens",
    "Lamellibrachia satsuma",
    "Lipotes vexillifer",
    "Lontra canadensis",
    "Loxodonta africana",
    "Macaca mulatta",
    "Monodelphis domestica",
    "Mus musculus",
    "Nematostella vectensis",
    "Notamacropus eugenii",
    "Ochotona princeps",
    "Octodon degus",
    "Ornithorhynchus anatinus",
    "Oryctolagus cuniculus",
    "Paramuricea clavata",
    "Pelodiscus sinensis",
    "Phascolarctos cinereus",
    "Physeter macrocephalus",
    "Podarcis muralis",
    "Salmo salar",
    "Sarcophilus harrisii",
    "Strongylocentrotus purpuratus",
    "Stylophora pistillata",
    "Sus scrofa",
    "Tursiops truncatus",
    "Ursus arctos",
    "Varanus komodoensis",
    "Vombatus ursinus",
    "Mirounga angustirostris",
    "Myotis lucifugus",
    "Oncorhynchus mykiss",
    "Orchesella cincta"
    ]

    # 2. Define your input and output filenames.
    INPUT_TSV_FILE = 'parsed_catsper_data.tsv'
    OUTPUT_CSV_FILE = 'catsper_length_summary_cleaned.csv'
    # ---------------------

    try:
        import pandas
    except ImportError:
        print("Error: The 'pandas' library is required. Please run: pip install pandas")
    else:
        create_cleaned_length_summary(INPUT_TSV_FILE, SPECIES_TO_INCLUDE, OUTPUT_CSV_FILE)