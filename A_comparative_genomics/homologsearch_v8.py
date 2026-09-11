import requests
import sys
import re
import time
import os
import xml.etree.ElementTree as ET

# --- Configuration ---
UNIPROT_API_URL = "https://rest.uniprot.org"
NCBI_API_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_API_KEY = "YOUR API KEY HERE"
DELAY_BETWEEN_REQUESTS = 0.11 # Adjust if no API key

CATSPER_DESCRIPTIVE_NAME_MAP = { # Maps primary ID to its descriptive search term part
    "CATSPER1": "cation channel sperm-associated protein 1",
    "CATSPER2": "cation channel sperm-associated protein 2",
    "CATSPER3": "cation channel sperm-associated protein 3",
    "CATSPER4": "cation channel sperm-associated protein 4",
    "CATSPERB": "cation channel sperm-associated protein subunit beta", # or "catsper channel auxiliary subunit beta"
    "CATSPERG": "cation channel sperm-associated protein subunit gamma",
    "CATSPERD": "cation channel sperm-associated protein subunit delta", # or "catsper channel auxiliary subunit delta"
    "CATSPERE": "cation channel sperm-associated protein subunit epsilon",
    "CATSPERZ": "cation channel sperm-associated protein subunit zeta",
    "EFCAB9": "EF-hand calcium-binding domain-containing protein 9",
    # TMEMs and SLC6A1 have less consistent long names and are usually found by their gene symbol.
}

CATSPER_TARGET_TERMS_BASE = [ # Original short names
    "CATSPER1", "CATSPER2", "CATSPER3", "CATSPER4", "CATSPERB",
    "CATSPERG", "CATSPERD", "CATSPERE", "CATSPERZ",
    "EFCAB9",
    "TMEM249", "TMEM262", "SLC6A1"
]
# Add descriptive names to the search list
CATSPER_TARGET_TERMS = list(CATSPER_TARGET_TERMS_BASE) # Start with a copy
for primary_id, desc_name in CATSPER_DESCRIPTIVE_NAME_MAP.items():
    if desc_name not in CATSPER_TARGET_TERMS: # Avoid duplicates if base already covers it
        CATSPER_TARGET_TERMS.append(desc_name)


# Ensure keys and values are uppercase for consistency
ALIAS_MAP_BASE = {
    "C2CD6": "EFCAB9",
    # Add mappings from descriptive names back to primary IDs
}
for primary_id, desc_name in CATSPER_DESCRIPTIVE_NAME_MAP.items():
    ALIAS_MAP_BASE[desc_name.upper()] = primary_id.upper()

ALIAS_MAP = {k.upper(): v.upper() for k, v in ALIAS_MAP_BASE.items()}

# PRIMARY_TARGETS are the canonical names we want in the final output
PRIMARY_TARGETS = sorted(list(set(ALIAS_MAP.get(t.upper(), t.upper()) for t in CATSPER_TARGET_TERMS_BASE)))


REQUEST_TIMEOUT = 60
REQUEST_SIZE = 30 # Number of results to fetch per query
DEBUG_PARSING = True

# --- Subunit Specificity Configuration (from previous version) ---
SUBUNIT_TYPE_KEYWORDS = {
    " 1": ["one", "alpha"], " 2": ["two", "beta"], " 3": ["three"], " 4": ["four"],
    "B": ["b", "beta"], "G": ["g", "gamma"], "D": ["d", "delta"],
    "E": ["e", "epsilon"], "Z": ["z", "zeta"],
}
PRIMARY_TARGET_TO_SPECIFIC_IDENTIFIER = {}
for pt in PRIMARY_TARGETS:
    if pt.startswith("CATSPER") and len(pt) > 7:
        identifier = pt[7:]
        if identifier in SUBUNIT_TYPE_KEYWORDS or identifier.isdigit():
            PRIMARY_TARGET_TO_SPECIFIC_IDENTIFIER[pt] = identifier

_sorted_primary_targets_for_regex = sorted(list(PRIMARY_TARGETS), key=len, reverse=True)
ANY_PRIMARY_TARGET_NAME_REGEX = re.compile(
    r"\b(" + "|".join(re.escape(t) for t in _sorted_primary_targets_for_regex) + r")(-like|-related\s+to|-type)?\b",
    re.IGNORECASE
)

all_word_kws = []
for type_kws_list in SUBUNIT_TYPE_KEYWORDS.values(): all_word_kws.extend(type_kws_list)
if all_word_kws:
    unique_word_kws = sorted(list(set(all_word_kws)), key=len, reverse=True)
    WORD_SUBUNIT_KEYWORD_REGEX_PATTERN = r"\b(" + "|".join(re.escape(kw) for kw in unique_word_kws) + r")\b"
    WORD_SUBUNIT_KEYWORD_REGEX = re.compile(WORD_SUBUNIT_KEYWORD_REGEX_PATTERN, re.IGNORECASE)
    if DEBUG_PARSING: print(f"DEBUG: Word Subunit Keyword Regex Pattern: {WORD_SUBUNIT_KEYWORD_REGEX_PATTERN}")
else:
    WORD_SUBUNIT_KEYWORD_REGEX = None

CONTEXTUAL_DIGIT_MARKERS = {}
POSSIBLE_CONTEXT_DIGITS = [" 1", " 2", " 3", " 4"]
DIGIT_CONTEXT_PREFIXES = ["protein", "subunit", "channel", "member", "type", "isoform", "component", "unit", "catsper", "sperm-associated protein"]
DIGIT_CONTEXT_SUFFIXES = ["like", "type"]
for d_str in POSSIBLE_CONTEXT_DIGITS:
    patterns = []
    for prefix in DIGIT_CONTEXT_PREFIXES:
        patterns.append(r"\b" + re.escape(prefix) + r"\s*[-_]?\s*" + re.escape(d_str) + r"(?![0-9a-zA-Z-])")
    for suffix in DIGIT_CONTEXT_SUFFIXES:
        patterns.append(r"(?<![0-9a-zA-Z-])" + re.escape(d_str) + r"-" + re.escape(suffix) + r"\b")
    # Add a pattern for just the digit if it's preceded by a space and followed by a space or -like (more general)
    patterns.append(r"(?<![a-zA-Z0-9-])" + re.escape(d_str) + r"(?:\s|\b|-like|-type)")


    if patterns:
        CONTEXTUAL_DIGIT_MARKERS[d_str] = re.compile("|".join(patterns), re.IGNORECASE)
        if DEBUG_PARSING:
            print(f"DEBUG: Contextual regex for digit marker '{d_str}': {CONTEXTUAL_DIGIT_MARKERS[d_str].pattern}")


# this is all the same
def execute_http_query(url: str, params: dict, term: str, db_name: str, method: str = "GET") -> str | None:
    try:
        if method.upper() == "POST": response = requests.post(url, data=params, timeout=REQUEST_TIMEOUT)
        else: response = requests.get(url, params=params, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
        return response.text
    except requests.exceptions.RequestException as e:
        print(f"\n    Error fetching {term} from {db_name}: {e}", file=sys.stderr)
    except requests.exceptions.HTTPError as e:
        if e.response.status_code != 404 or DEBUG_PARSING:
             print(f"\n    HTTP Error fetching {term} from {db_name}: {e}. Status: {e.response.status_code}", file=sys.stderr)
    return None

def _extract_header_from_fasta_entry(fasta_entry: str) -> str:
    return fasta_entry.split('\n', 1)[0] if fasta_entry else ""

def _extract_sequence_from_fasta_entry(fasta_entry: str) -> str:
    if not fasta_entry: return ""
    lines = fasta_entry.strip().split('\n')
    return "".join(s.strip().upper() for s in lines[1:]) if len(lines) > 1 else ""

def _split_multifasta(fasta_string: str) -> list[str]:
    if not fasta_string or not fasta_string.strip().startswith('>'): return []
    entries = fasta_string.replace('\r\n', '\n').replace('\r', '\n').strip().split('\n>')
    return [entries[0]] + ['>' + part for part in entries[1:]] if entries else []

def construct_uniprotkb_query(taxid: int, term: str) -> str: return f'("{term}" AND taxonomy_id:{taxid})'
def construct_uniparc_query(taxid: int, term: str) -> str: return f'("{term}" AND organism_id:{taxid})' # UniParc often uses organism_id

def fetch_raw_candidates_from_uniprot(search_endpoint: str, params: dict, term: str, db_name: str) -> list[str]:
    raw_fasta_data = execute_http_query(f"{UNIPROT_API_URL}{search_endpoint}", params, term, db_name)
    return _split_multifasta(raw_fasta_data) if raw_fasta_data else []

def fetch_raw_candidates_from_ncbi(taxid: int, term: str, request_size: int) -> list[str]:
    # Broader search for NCBI: include the term in [All Fields] as well
    # but prioritize [Protein Name] or [Gene Name] if possible by how NCBI ranks results.
    # The query structure here sends one query with ORs.
    search_term_for_ncbi = f'("{term}"[Protein Name] OR "{term}"[Gene Name] OR "{term}"[Title])'
    # If term is very generic, [All Fields] can be too noisy. Title is a good compromise.
    # If term itself is a longer descriptive phrase, it might directly match the title.
    
    esearch_params = { "db": "protein",
                       "term": f'{search_term_for_ncbi} AND txid{taxid}[Organism:exp]',
                       "retmax": str(request_size), "retmode": "xml", "usehistory": "n"}
    if NCBI_API_KEY: esearch_params["api_key"] = NCBI_API_KEY
    esearch_xml = execute_http_query(f"{NCBI_API_URL}esearch.fcgi", esearch_params, term, "NCBI ESearch")
    time.sleep(DELAY_BETWEEN_REQUESTS)
    if not esearch_xml: return []
    try: uids = [id_tag.text for id_tag in ET.fromstring(esearch_xml).findall(".//Id")]
    except ET.ParseError as e: print(f"\n    Error parsing NCBI ESearch XML for '{term}': {e}", file=sys.stderr); return []
    if not uids: return []
    
    efetch_params = {"db": "protein", "id": ",".join(uids), "rettype": "fasta", "retmode": "text"}
    if NCBI_API_KEY: efetch_params["api_key"] = NCBI_API_KEY
    raw_fasta = execute_http_query(f"{NCBI_API_URL}efetch.fcgi", efetch_params, term, "NCBI EFetch")
    return _split_multifasta(raw_fasta) if raw_fasta else []

def get_taxid_from_user() -> int:
    while True:
        try:
            taxid = int(input("Enter the NCBI Taxonomy ID: "))
            if taxid > 0: return taxid
            print("TaxID must be a positive integer.")
        except ValueError: print("Invalid input. Please enter an integer.")
        except EOFError: print("\nInput cancelled. Exiting.", file=sys.stderr); sys.exit(1)

def get_allowed_keywords_for_target_type(target_specific_id: str | None) -> set[str]:
    if not target_specific_id: return set()
    allowed_kws = set()
    allowed_kws.update(kw.lower() for kw in SUBUNIT_TYPE_KEYWORDS.get(target_specific_id, []))
    if target_specific_id.isdigit(): # If target_specific_id is "1", "1" is allowed.
        allowed_kws.add(target_specific_id)
    return allowed_kws

def is_header_acceptable_for_target(header: str, intended_primary_id: str) -> bool:
    header_lower = header.lower()

    # Check 1: Explicit mention of a *different* primary target name
    for match_obj in ANY_PRIMARY_TARGET_NAME_REGEX.finditer(header):
        mentioned_target_full = match_obj.group(1).upper()
        mentioned_target_primary = ALIAS_MAP.get(mentioned_target_full, mentioned_target_full)
        if mentioned_target_primary != intended_primary_id and mentioned_target_primary in PRIMARY_TARGETS:
            if DEBUG_PARSING: 
                print(f"        REJECT (DiffPrim): Mentions primary '{mentioned_target_primary}' (matched '{match_obj.group(0)}') for intended '{intended_primary_id}'. H: {header[:90]}")
            return False

    intended_subunit_specific_id = PRIMARY_TARGET_TO_SPECIFIC_IDENTIFIER.get(intended_primary_id)
    allowed_kws_for_intended_type = get_allowed_keywords_for_target_type(intended_subunit_specific_id)

    # Check 2: Mismatch with WORD-BASED subunit type keywords
    if WORD_SUBUNIT_KEYWORD_REGEX:
        for match_obj in WORD_SUBUNIT_KEYWORD_REGEX.finditer(header_lower):
            found_word_keyword = match_obj.group(1).lower() # group(1) is the captured keyword
            is_a_defined_word_keyword = any(found_word_keyword in type_list for type_list in SUBUNIT_TYPE_KEYWORDS.values()) # Check against original SUBUNIT_TYPE_KEYWORDS lists
            if is_a_defined_word_keyword and found_word_keyword not in allowed_kws_for_intended_type:
                if DEBUG_PARSING:
                    print(f"        REJECT (WordSubunitKw): Mentions word_kw '{found_word_keyword}' (matched '{match_obj.group(0)}') not in allowed {allowed_kws_for_intended_type} for type '{intended_subunit_specific_id}' of '{intended_primary_id}'. H: {header[:90]}")
                return False

    # Check 3: Mismatch with DIGIT identifiers in specific contexts
    for digit_marker_str, contextual_regex_for_this_digit_marker in CONTEXTUAL_DIGIT_MARKERS.items():
        # Extract the core digit/identifier from the marker string (e.g., ' 1' -> '1')
        core_digit_marker = digit_marker_str.strip()
        is_this_digit_marker_allowed = core_digit_marker in allowed_kws_for_intended_type
        
        match_in_context = contextual_regex_for_this_digit_marker.search(header_lower)
        if match_in_context:
            if not is_this_digit_marker_allowed: # Found a contextual digit that is NOT allowed for current target
                if DEBUG_PARSING:
                    print(
                        f"        REJECT (ContextualDigit): Header contains contextual unwanted digit '{digit_marker_str}' "
                        f"(matched '{match_in_context.group(0)}' at {match_in_context.span()}) "
                        f"for intended '{intended_primary_id}' (type: '{intended_subunit_specific_id}', allowed kws: {allowed_kws_for_intended_type}). H: {header[:90]}"
                    )
                return False
            elif DEBUG_PARSING: # Found a contextual digit that IS allowed - this is fine.
                 print(
                        f"        INFO (ContextualDigitOK): Header contains contextual allowed digit '{digit_marker_str}' "
                        f"(matched '{match_in_context.group(0)}' at {match_in_context.span()}) "
                        f"for intended '{intended_primary_id}'. This specific match is OK.")
    return True

# --- Main function and subsequent phases (Phase 1, 2, 3, Writing) remain the same ---
def main():
    print("--- CatSper Sequence Retrieval (v8: Broader Search Terms) ---")
    if NCBI_API_KEY: print("Using NCBI API Key.")
    else: print("No NCBI API Key. Rate limits: max ~3 requests/sec. Delays enforced.")
    if not NCBI_API_KEY and DELAY_BETWEEN_REQUESTS < 0.33:
        print(f"Warning: DELAY_BETWEEN_REQUESTS ({DELAY_BETWEEN_REQUESTS}s) is too short for NCBI without API key. Adjusting to 0.34s.")
        globals()['DELAY_BETWEEN_REQUESTS'] = 0.34

    print(f"Will fetch up to {REQUEST_SIZE} entries from each DB per term.")
    print(f"Search terms to be used: {', '.join(CATSPER_TARGET_TERMS)}")


    taxid = get_taxid_from_user()
    output_dir_name = "results"
    if not os.path.exists(output_dir_name):
        os.makedirs(output_dir_name); print(f"Created output directory: {output_dir_name}")
    base_filename = f"catsper_v8_{taxid}.fasta"
    output_file_path = os.path.join(output_dir_name, base_filename)

    print(f"\nAttempting retrieval for TaxID: {taxid}")
    # print(f"Target Search Terms: {', '.join(CATSPER_TARGET_TERMS)}") # Already printed above
    print(f"Primary Targets (for final output): {', '.join(PRIMARY_TARGETS)}")
    print(f"Output will be saved to: {output_file_path}")

    # --- Phase 1: Fetch ALL potential candidates ---
    all_potential_candidates_by_primary_id = {pid: [] for pid in PRIMARY_TARGETS}
    # Keep track of which (term_to_search, taxid) combinations we've already queried to avoid redundant API calls
    # if multiple CATSPER_TARGET_TERMS map to the same primary_id and would generate identical API queries.
    # However, different terms (e.g. "CATSPER1" vs "cation channel... protein 1") WILL generate different queries.
    # The `processed_search_terms` in the previous version was to skip an alias if the primary id had been searched.
    # Here, we want to search for ALL terms in CATSPER_TARGET_TERMS.
    # The grouping happens when we store results into `all_potential_candidates_by_primary_id`.

    print("\n--- Phase 1: Fetching all potential candidates ---")
    for term_to_search in CATSPER_TARGET_TERMS: # Iterate through the expanded list
        primary_intended_id = ALIAS_MAP.get(term_to_search.upper(), term_to_search.upper())
        
        # Ensure primary_intended_id is one of our known canonical targets
        if primary_intended_id not in PRIMARY_TARGETS:
            if DEBUG_PARSING:
                print(f"  Skipping search term '{term_to_search}' as its primary ID '{primary_intended_id}' is not in the main PRIMARY_TARGETS list.")
            continue

        print(f"  Fetching for search term '{term_to_search}' (maps to primary: {primary_intended_id})")
        current_term_raw_entries = []
        # UniProtKB
        print(f"    Querying UniProtKB... ", end="")
        kb_params = {'query': construct_uniprotkb_query(taxid, term_to_search), 'format': 'fasta', 'size': REQUEST_SIZE, 'includeIsoforms': 'true'}
        kb_raw = fetch_raw_candidates_from_uniprot("/uniprotkb/search", kb_params, term_to_search, "UniProtKB")
        current_term_raw_entries.extend([(e, "UniProtKB") for e in kb_raw]); print(f"Got {len(kb_raw)}.")
        time.sleep(DELAY_BETWEEN_REQUESTS)
        # UniParc
        print(f"    Querying UniParc... ", end="")
        parc_params = {'query': construct_uniparc_query(taxid, term_to_search), 'format': 'fasta', 'size': REQUEST_SIZE}
        parc_raw = fetch_raw_candidates_from_uniprot("/uniparc/search", parc_params, term_to_search, "UniParc")
        current_term_raw_entries.extend([(e, "UniParc") for e in parc_raw]); print(f"Got {len(parc_raw)}.")
        time.sleep(DELAY_BETWEEN_REQUESTS)
        # NCBI
        print(f"    Querying NCBI... ", end="")
        ncbi_raw = fetch_raw_candidates_from_ncbi(taxid, term_to_search, REQUEST_SIZE)
        current_term_raw_entries.extend([(e, "NCBI") for e in ncbi_raw]); print(f"Got {len(ncbi_raw)}.")
        time.sleep(DELAY_BETWEEN_REQUESTS) # Delay after the group of fetches for this term_to_search

        for fasta_entry_str, source_db_str in current_term_raw_entries:
            seq_str = _extract_sequence_from_fasta_entry(fasta_entry_str)
            if seq_str:
                # Store under the canonical primary_intended_id
                all_potential_candidates_by_primary_id[primary_intended_id].append(
                    (fasta_entry_str, seq_str, len(seq_str), source_db_str, term_to_search)
                )

    # --- Phase 2: Filter candidates and select initial best ---
    print("\n--- Phase 2: Initial selection of longest, acceptable candidate per primary target ---")
    initial_best_selection = {}
    valid_candidates_for_dedup = {pid: [] for pid in PRIMARY_TARGETS}

    for primary_id in PRIMARY_TARGETS: # Iterate through canonical primary target names
        # Consolidate and unique-ify candidates for this primary_id from potentially multiple search terms
        unique_raw_candidates_for_pid_by_header = {} # Use header to get unique entries before filtering
        for cand_tuple in all_potential_candidates_by_primary_id.get(primary_id, []):
            header = _extract_header_from_fasta_entry(cand_tuple[0])
            if header not in unique_raw_candidates_for_pid_by_header:
                 unique_raw_candidates_for_pid_by_header[header] = cand_tuple
        
        consolidated_candidates = list(unique_raw_candidates_for_pid_by_header.values())

        if not consolidated_candidates:
            if DEBUG_PARSING: print(f"  No raw candidates found for {primary_id} from any search term.")
            continue
        
        if DEBUG_PARSING: print(f"  Filtering and selecting for {primary_id} from {len(consolidated_candidates)} unique raw candidates:")
        
        current_pid_acceptable_options = []
        for fasta_entry, seq_str, seq_len, source_db, original_search_term in consolidated_candidates:
            header = _extract_header_from_fasta_entry(fasta_entry) # Already got it, but for clarity
            if is_header_acceptable_for_target(header, primary_id):
                current_pid_acceptable_options.append((fasta_entry, seq_str, seq_len, source_db, original_search_term))
                # Add to valid_candidates_for_dedup which is {primary_id: [list_of_valid_tuples]}
                valid_candidates_for_dedup[primary_id].append((fasta_entry, seq_str, seq_len, source_db, original_search_term))
                if DEBUG_PARSING: print(f"    [+] ACCEPTED for {primary_id} (orig search: '{original_search_term}'). Len: {seq_len}, Src: {source_db}. H: {header[:80]}...")
            
        if current_pid_acceptable_options:
            current_pid_acceptable_options.sort(key=lambda x: x[2], reverse=True)
            best_opt = current_pid_acceptable_options[0]
            initial_best_selection[primary_id] = (best_opt[0], best_opt[1], best_opt[3], best_opt[4], best_opt[2])
            print(f"  Initial best for {primary_id}: Length {best_opt[2]} from {best_opt[3]} (orig search: '{best_opt[4]}').")
        else:
            print(f"  No acceptable candidates found for {primary_id} after filtering.")

    # --- Step 3: De-duplication and final selection ---
    print("\n--- Phase 3: De-duplication and Final Selection ---")
    final_selected_sequences = {} 
    used_sequence_strings = set()

    for primary_id in PRIMARY_TARGETS: 
        if primary_id not in initial_best_selection:
            if DEBUG_PARSING: print(f"  Skipping de-dup for {primary_id}: No initial selection.")
            continue

        current_best_fasta, current_best_seq, current_best_src, _, current_best_len = initial_best_selection[primary_id]

        if current_best_seq not in used_sequence_strings:
            final_selected_sequences[primary_id] = (current_best_fasta, current_best_src)
            used_sequence_strings.add(current_best_seq)
            print(f"  Selected for {primary_id}: Original best (Length {current_best_len}, Source {current_best_src}). Sequence is unique.")
        else:
            print(f"  Conflict for {primary_id}: Initial best sequence (Length {current_best_len}) already used. Looking for alternative...")
            # Get unique alternatives from valid_candidates_for_dedup[primary_id]
            # Ensure alternatives are unique by sequence before sorting by length
            unique_alt_sequences = {}
            for alt_tuple in valid_candidates_for_dedup.get(primary_id, []):
                alt_seq = alt_tuple[1]
                if alt_seq not in unique_alt_sequences:
                    unique_alt_sequences[alt_seq] = alt_tuple
            
            alternatives_sorted = sorted(list(unique_alt_sequences.values()), key=lambda x: x[2], reverse=True)
            
            found_alternative = False
            for alt_fasta, alt_seq_str, alt_len, alt_source_db, _ in alternatives_sorted:
                if alt_seq_str not in used_sequence_strings:
                    final_selected_sequences[primary_id] = (alt_fasta, alt_source_db)
                    used_sequence_strings.add(alt_seq_str)
                    print(f"    Found alternative for {primary_id}: Length {alt_len}, Source {alt_source_db}. Sequence is unique.")
                    found_alternative = True
                    break
            if not found_alternative:
                print(f"    No unique alternative sequence found for {primary_id}. It will be excluded.")
    
    # --- Writing Results ---
    print(f"\nCollected {len(final_selected_sequences)} final unique sequences.")
    if not final_selected_sequences:
        print("\n--- No sequences selected after all phases. ---")
    
    if final_selected_sequences:
        try:
            with open(output_file_path, 'w') as f:
                f.write(f"# Sequences retrieved for TaxID: {taxid}\n")
                f.write(f"# Selection: Longest, header-acceptable, unique sequence per primary target.\n")
                f.write(f"# Fetched up to {REQUEST_SIZE} candidates per DB per search term initially.\n")
                f.write(f"# Format: >PrimaryTargetID | Source Database=SOURCE | Original Header...\n\n")
                for pid_key in sorted(final_selected_sequences.keys()):
                    fasta_entry, source = final_selected_sequences[pid_key]
                    original_header_line = _extract_header_from_fasta_entry(fasta_entry)
                    sequence_part = _extract_sequence_from_fasta_entry(fasta_entry)
                    
                    header_content = original_header_line.lstrip('>').strip()
                    new_header = f">{pid_key} | Source Database={source} | {header_content}"
                    
                    f.write(new_header + "\n")
                    f.write(sequence_part + "\n") 
                    f.write("\n")
            print(f"Successfully wrote {len(final_selected_sequences)} unique sequences to: {output_file_path}")
        except IOError as e:
            print(f"Error writing output file {output_file_path}: {e}", file=sys.stderr)
            sys.exit(1)

    configured_primary_targets_set = set(PRIMARY_TARGETS)
    found_primary_ids_in_final = set(final_selected_sequences.keys())
    missing_primary_in_final = configured_primary_targets_set - found_primary_ids_in_final

    if missing_primary_in_final:
        print(f"Warning: Could not select unique sequences for: {', '.join(sorted(list(missing_primary_in_final)))}")
    elif final_selected_sequences:
        print("Found representative unique sequences for all expected primary target terms.")

    print("\n--- Process Completed ---")

if __name__ == "__main__":
    main()





    