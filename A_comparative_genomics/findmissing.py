import os
from pathlib import Path
from typing import List

def find_fasta_without_phrase(
    folder_path: str, 
    search_phrase: str, 
    case_sensitive: bool = False
) -> List[str]:
    """
    Scans a folder for .fasta files and reports which ones DO NOT contain a specific phrase.

    Args:
        folder_path (str): The path to the folder containing your .fasta files.
        search_phrase (str): The phrase to search for within each file.
        case_sensitive (bool): Set to True for a case-sensitive search. 
                               Defaults to False (recommended).

    Returns:
        List[str]: A list of filenames that are missing the search_phrase. 
                   Returns an empty list if all files contain the phrase.
    """
    # --- 1. Validate the folder path ---
    folder = Path(folder_path)
    if not folder.is_dir():
        print(f"Error: The folder '{folder_path}' does not exist.")
        return []

    # --- 2. Find all relevant FASTA files ---
    # You can add more extensions here if needed, like '*.fna' or '*.ffn'
    fasta_extensions = ['*.fasta', '*.fa']
    fasta_files = []
    for ext in fasta_extensions:
        fasta_files.extend(folder.glob(ext))

    if not fasta_files:
        print(f"Warning: No .fasta or .fa files found in '{folder_path}'.")
        return []

    # --- 3. Check each file for the phrase ---
    files_missing_phrase = []
    
    # Adjust search phrase for case-insensitivity if needed
    phrase_to_find = search_phrase if case_sensitive else search_phrase.lower()

    for file_path in fasta_files:
        print(f"Checking {file_path.name}...")
        found = False
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                # Read file line-by-line for memory efficiency
                for line in f:
                    line_to_check = line if case_sensitive else line.lower()
                    if phrase_to_find in line_to_check:
                        found = True
                        break # Phrase found, no need to read the rest of the file
            
            if not found:
                files_missing_phrase.append(file_path.name)
        
        except Exception as e:
            print(f"  -> Could not read file {file_path.name}. Error: {e}")
            # You might want to add files that can't be read to a separate list
    
    return files_missing_phrase


if __name__ == "__main__":
    # SET FOLDER AND PHRASE HERE
    folder_with_fasta = "./results/"  # "." means the current directory where the script is running
    phrase_to_look_for = ">TMEM262" # Example phrase

    print(f"Searching for files in '{os.path.abspath(folder_with_fasta)}' that DON'T contain the phrase: '{phrase_to_look_for}'\n")

    # The search is case-insensitive by default. To make it sensitive, add `case_sensitive=True`
    missing_list = find_fasta_without_phrase(folder_with_fasta, phrase_to_look_for)

    # 3. PRINT THE RESULTS
    print("\n--- Scan Complete ---")
    if not missing_list:
        print("All FASTA files contain the specified phrase.")
    else:
        print(f"Found {len(missing_list)} file(s) that DO NOT contain the phrase:")
        for filename in missing_list:
            print(f"  - {filename}")