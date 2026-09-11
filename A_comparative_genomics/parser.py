import os
import re
import csv
from ete3 import NCBITaxa # Import the NCBI Taxonomy tool

# --- ETE3 Initialization ---
# This creates an interface to the NCBI Taxonomy database.
# IMPORTANT: The first time you run this script, ete3 will download the
# taxonomy database. This may take a few minutes and requires an internet connection.
# Subsequent runs will be fast and can work offline.
print("Initializing NCBI Taxonomy database (may download on first run)...")
ncbi = NCBITaxa()
print("Database ready.")

# A cache to store TaxID -> Species Name mappings to avoid redundant lookups
taxid_cache = {}

def get_species_name(taxid, header):
    """
    Finds the species name. First, it tries parsing the header for an 'OS=' tag.
    If that fails, it uses the taxid to query the NCBI database.
    """
    # 1. Try the fast way: parse the header
    species_pattern = re.compile(r'OS=([^=]+?)(?:\sOX=|\sGN=|\sPE=|\sSV=|$)')
    match = species_pattern.search(header)
    if match:
        return match.group(1).strip()

    # 2. If parsing fails, use the fallback: NCBI lookup via TaxID
    # Check our cache first to see if we've already looked this up
    if taxid in taxid_cache:
        return taxid_cache[taxid]
    
    # If not in cache, query NCBI, then store the result
    try:
        taxid_int = int(taxid)
        # get_taxid_translator wants a list of IDs
        translation = ncbi.get_taxid_translator([taxid_int])
        if translation:
            species_name = translation[taxid_int]
            taxid_cache[taxid] = species_name # Store in cache
            return species_name
    except Exception as e:
        print(f"    - Warning: Could not look up TaxID {taxid}. Error: {e}")
        
    # If all else fails, return a default value
    taxid_cache[taxid] = 'Unknown Species'
    return 'Unknown Species'


def parse_fasta_files(input_folder, output_file):
    """
    Parses FASTA files, intelligently finds the species name, and creates a TSV.
    """
    filename_pattern = re.compile(r'catsper_v8_(\d+)\.fasta')
    if not os.path.isdir(input_folder):
        print(f"Error: Input folder '{input_folder}' not found.")
        return

    print(f"Starting to process files in '{input_folder}'...")

    with open(output_file, 'w', newline='') as f_out:
        writer = csv.writer(f_out, delimiter='\t')
        writer.writerow([
            "TaxID", "SpeciesName", "GeneName", "SourceDatabase",
            "OriginalHeader", "SequenceLength", "Sequence"
        ])

        for filename in sorted(os.listdir(input_folder)):
            match = filename_pattern.match(filename)
            if not match:
                continue

            taxid = match.group(1)
            file_path = os.path.join(input_folder, filename)
            print(f"  -> Processing file: {filename} (TaxID: {taxid})")

            with open(file_path, 'r') as f_in:
                header = None
                sequence_parts = []
                for line in f_in:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    if line.startswith('>'):
                        if header:
                            full_sequence = "".join(sequence_parts)
                            process_record(writer, taxid, header, full_sequence)
                        header = line
                        sequence_parts = []
                    else:
                        if header:
                            sequence_parts.append(line)
                if header:
                    full_sequence = "".join(sequence_parts)
                    process_record(writer, taxid, header, full_sequence)

    print(f"\nProcessing complete. Data saved to '{output_file}'.")


def process_record(writer, taxid, full_header, sequence):
    """Parses a single record and writes it to the TSV."""
    try:
        header_content = full_header[1:]
        
        # Use our new, smart function to get the species name
        species_name = get_species_name(taxid, header_content)

        parts = [p.strip() for p in header_content.split('|')]
        gene_name = parts[0]
        source_db = parts[1].split('=', 1)[1] if len(parts) > 1 and '=' in parts[1] else 'N/A'
        original_header = parts[2].strip() if len(parts) > 2 else 'N/A'
        
        writer.writerow([
            taxid, species_name, gene_name, source_db,
            original_header, len(sequence), sequence
        ])
    except Exception as e:
        print(f"    - Error processing record with header '{full_header}': {e}")


if __name__ == '__main__':
    INPUT_FOLDER = 'results' 
    OUTPUT_FILE = 'parsed_catsper_data.tsv'
    parse_fasta_files(INPUT_FOLDER, OUTPUT_FILE)