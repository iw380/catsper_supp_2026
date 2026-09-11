import os
import subprocess

# --- Configuration ---
INPUT_FASTA = "catsper1_n_terminals.fasta"
OUTPUT_ALIGNMENT = "aligned_nterms_filtered.fasta"
FILTERED_FASTA = "filtered_sequences.fasta"

# Species you want to keep (exactly match or partial match)
SPECIES_TO_KEEP = [
    "Acropora_digitifera",
    "Acropora_millepora",
    "Alligator_mississippiensis",
    "Apostichopus_japonicus",
    "Bos_taurus",
    "Canis_lupus_familiaris",
    "Capra_hircus",
    "Capricornis_sumatraensis",
    "Cavia_porcellus",
    "Ciona_intestinalis",
    "Elephantulus_edwardii",
    "Equus_caballus",
    "Fukomys_damarensis",
    "Heterocephalus_glaber",
    "Homo_sapiens",
    "Lamellibrachia_satsuma",
    "Lingula_anatina",
    "Lipotes_vexillifer",
    "Lontra_canadensis",
    "Loxodonta_africana",
    "Macaca_mulatta",
    "Mirounga_angustirostris",
    "Monodelphis_domestica",
    "Mus_musculus",
    "Myotis_lucifugus",
    "Nematostella_vectensis",
    "Notamacropus_eugenii",
    "Ochotona_princeps",
    "Octodon_degus",
    "Oncorhynchus_mykiss",
    "Orchesella_cincta",
    "Ornithorhynchus_anatinus",
    "Oryctolagus_cuniculus",
    "Paramuricea_clavata",
    "Pelodiscus_sinensis",
    "Phascolarctos_cinereus",
    "Physeter_macrocephalus",
    "Podarcis_muralis",
    "Salmo_salar",
    "Sarcophilus_harrisii",
    "Strongylocentrotus_purpuratus",
    "Stylophora_pistillata",
    "Sus_scrofa",
    "Tursiops_truncatus",
    "Ursus_arctos",
    "Varanus_komodoensis",
    "Vombatus_ursinus"
]

# ---------------------

def filter_fasta_by_species(input_file, output_file, allowed_species):
    """
    Reads a FASTA file and writes a new FASTA containing only entries
    whose headers contain one of the allowed species names.
    """
    if not os.path.exists(input_file):
        print(f"--- ERROR: Input file not found at '{input_file}' ---")
        return None

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        keep = False
        kept_count = 0

        for line in infile:
            if line.startswith(">"):
                # Check if any allowed species name appears in the header
                keep = any(species in line for species in allowed_species)
                if keep:
                    outfile.write(line)
                    kept_count += 1
            else:
                if keep:
                    outfile.write(line)

    if kept_count == 0:
        print("--- WARNING: No sequences matched the allowed species list! ---")
    else:
        print(f"Filtered FASTA written to '{output_file}' with {kept_count} species included.")

    return output_file


def run_muscle_alignment(input_file, output_file):
    """
    Runs a MUSCLE multiple sequence alignment on an input FASTA file
    and saves the result to an output file.
    """
    if not os.path.exists(input_file):
        print(f"--- ERROR: Input file not found at '{input_file}' ---")
        return

    print(f"Starting multiple sequence alignment for: '{input_file}'")

    try:
        subprocess.run(
            ["muscle", "-align", input_file, "-output", output_file],
            check=True,
            capture_output=True,
            text=True
        )
        print(f"--- Alignment successful! ---")
        print(f"Aligned sequences have been saved to: '{output_file}'")

    except FileNotFoundError:
        print("\n--- MUSCLE ERROR ---")
        print("The command 'muscle' was not found. Ensure MUSCLE is installed and in your PATH.")
    except subprocess.CalledProcessError as e:
        print("\n--- MUSCLE FAILED ---")
        print("Error details from MUSCLE:")
        print(e.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


# --- Main Execution ---
if __name__ == "__main__":
    filtered = filter_fasta_by_species(INPUT_FASTA, FILTERED_FASTA, SPECIES_TO_KEEP)
    if filtered:
        run_muscle_alignment(filtered, OUTPUT_ALIGNMENT)
