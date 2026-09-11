# Supplementary Files

Supplementary files for the 2026 Catsper paper. 

## A. Comparative genomics & sequence alignment
- `homologsearch_v8.py` — retrieves sequences via UniProtKB/UniParc/NCBI APIs
- `results/` — per-taxon FASTA output of `homologsearch_v8.py`
- `parser.py` — parses `results/` into `parsed_catsper_data.tsv` (uses `ete3` for species-name resolution, backed by a local NCBI taxonomy dump not included here)
- `parsed_catsper_data.tsv` — master retrieval table
- `lister.py` + `species_list.txt` — unique species list (input to NCBI Common Tree)
- `findmissing.py` — QC check for missing gene coverage
- `csv_allsubunits.py` + `catsper_length_summary_cleaned.csv` — IQR-outlier-filtered per-subunit length summary (supports isoform/species curation)

## B. Definition of the CatSper1 N-Terminal Region
- `getnterminalranges.py` — MUSCLE-aligns full-length CATSPER1 sequences, finds first ≥10-column block with >80% occupancy → `catsper1_nterminal_ranges.tsv`
- `catsper1_nterminal_ranges.tsv`
- `extract_n_terminals.py` — slices N-terminal sequences using the same conserved-block rule into `catsper1_n_terminals.fasta`
- `catsper1_n_terminals.fasta`
- `get_alignment.py` — filters `catsper1_n_terminals.fasta` down to a list of 46 high-quality sequences into `filtered_sequences.fasta`, then runs `muscle -align` (MUSCLE v5.2) on that filtered set with default parameters into `aligned_nterms_filtered.fasta`
- `filtered_sequences.fasta`, `aligned_nterms_filtered.fasta`

## C. Quantification of N-Terminal Histidine Content
- `create_csv.py` — computes His count / N-terminal length (%) per species
- `master_dataset_primary.csv` — output, incl. `N_Terminal_Histidine_Percent`

## D. Physicochemical Properties (pI, LLPS)
- `manual_pi_values.txt` — per-species isoelectric points
- `profile_data/catsper[1-4]_fuzdrop.tsv` — raw FuzDrop pDP/Sbind profiles for the four pore-forming subunits (mouse reference)
- `sslp_plots.py` — generates the pDP/Sbind profile figures

## E. Catspermasome Structural Modeling & RSA/SASA Analysis
- `pdb_manifest.txt` — maps PDB filenames to species
- `surface_hist.py` — RSA/SASA calculation (RSA_THRESHOLD=25.0, Tien et al. 2013 MaxSASA table, Bio.PDB `SASA.ShrakeRupley`)
- `create_csv.py` — integrates the same RSA/SASA logic into the master table
- `master_dataset_surface_histidines.csv` — final per-species surface-histidine dataset

## F. NARDINI+ 2 Analysis
The Run_NARDINI+ directory extracts 90 feature z-score vectors given a fasta file of IDRs
#################### STEP 1: RUN NARDINI+ AND EXTRACT Z-SCORE VECTORS ################################

1. Run remove_bad_idrs.ipynb
   - This python script removes IDRs that are shorter than 30 and / or contain a character that isn't one of the twenty amino acids
   - The cutoff of 30 is used to avoid erroneous interpretations from short sequences when calculating features based on fractions of a given residue or residue type
   - The output of this script is a cleaned fasta file and a file with just one IDR sequence per line
   - The latter is used as input in NARDINI
   
2. Run submit2_slurm.bsub
   - This submits NARDINI jobs to a cluster using the python script nardini_calc_set_b1_random_seed.py
   - Individual NARDINI job files per sequence are output in the zscores directory
   - How you submit NARDINI jobs will be dependent on your cluster set up 
   - Even without a cluster you can just run the nardini_calc_set_b1_random_seed.py script where the output of interest is the zm_list
   
3. Run analyze_sequence_gramamr_from_FASTA.ipynb
   - This script extracts the full 90 feature z-score vectors
   - For compositional features we use the human IDRome as the prior which was previously extracted in Ruff et al. Cell 2026 Molecular grammars of predicted intrinsically disordered regions that span the human proteome
   - The main cells to run are Part 0, Part 1, 1.1, 1.2, 2.1, 2.2, 3.1, 3.2, 3.3, 4.1, and 4.2
   - 4.2 outputs the excel file of the full 90 feature z-score vectors per IDR
   
#################### STEP 2: ANALZYE NARDINI+ DATA ################################
The Analyze_NARDINI+ directory analyses NARDINI+ data and plots figures 2a-c

1. Run analyze_features_different_between_fertilization_temperatures.ipynb
   - This script determines which grammar features have different distributions for Cold, Intermediate, and Warm fertilizing species (Fig 2c)
   - Plots the distributions of selected grammar features that have different distributions (Fig 2b)
   - For the selected grammar features plots the z-scores per species (Fig 2a)
   - This script also outputs excel files that can used to plot Figures 2a, b, c
