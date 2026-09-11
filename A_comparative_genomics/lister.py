import pandas as pd

def extract_unique_species(input_tsv, output_txt):
    """
    Reads a TSV file, extracts all unique species names, sorts them,
    and writes them to a text file, one per line.

    Args:
        input_tsv (str): Path to the input TSV file from the parsing script.
        output_txt (str): Path to the output text file for the species list.
    """
    print(f"Reading data from '{input_tsv}'...")
    try:
        # Load the data using pandas
        df = pd.read_csv(input_tsv, sep='\t')
    except FileNotFoundError:
        print(f"--- ERROR ---")
        print(f"The input file '{input_tsv}' was not found.")
        print("Please make sure you have run 'parse_fasta.py' to generate it.")
        return

    # Check if the 'SpeciesName' column exists
    if 'SpeciesName' not in df.columns:
        print(f"--- ERROR ---")
        print(f"The column 'SpeciesName' was not found in '{input_tsv}'.")
        print("Please re-run the latest version of 'parse_fasta.py' to include it.")
        return

    # Get all unique species names from the column
    # The .unique() method returns an array of unique values
    unique_species = df['SpeciesName'].unique()

    # --- Clean up the list ---
    # Convert to a list, filter out any non-string or placeholder values
    placeholders = {'n/a', 'unknown species'}
    cleaned_list = [
        name for name in unique_species 
        if isinstance(name, str) and name.lower() not in placeholders
    ]

    # Sort the list alphabetically
    sorted_list = sorted(cleaned_list)

    # --- Write the list to the output file ---
    print(f"Found {len(sorted_list)} unique species.")
    with open(output_txt, 'w') as f_out:
        for species_name in sorted_list:
            f_out.write(species_name + '\n')

    print(f"Successfully saved the list to '{output_txt}'.")


if __name__ == '__main__':
    # --- Configuration ---
    INPUT_FILE = 'parsed_catsper_data.tsv'
    OUTPUT_FILE = 'species_list.txt'

    # Check for required libraries
    try:
        import pandas
    except ImportError:
        print("The 'pandas' library is required but not installed.")
        print("Please install it by running: pip install pandas")
    else:
        # --- Run the extraction ---
        extract_unique_species(INPUT_FILE, OUTPUT_FILE)