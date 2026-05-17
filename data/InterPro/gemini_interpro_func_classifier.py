import google.generativeai as genai
import json
import os
import time
import logging
import re # Import regular expressions for parsing
import csv
import datetime

# --- Basic Logging Configuration ---
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("interpro_classification.log"), # Log to a file
        logging.StreamHandler() # Also log to console
    ]
)

# --- Configuration ---
API_KEY = os.getenv("GOOGLE_API_KEY")
if not API_KEY:
    logging.error("GOOGLE_API_KEY environment variable not set.")
    raise ValueError("GOOGLE_API_KEY environment variable not set.")

# --->> !! CHANGE THESE !! <<---
JSON_FILE_PATH = "/Users/yangyxt/Downloads/Interpro_entry_mapping.json"  # Path to your input JSON file
OUTPUT_FILE_PATH = "classified_functional_entries.json" # Path to save results
# -----------------------------

BATCH_SIZE = 40
MODEL_NAME = "gemini-2.0-flash" # Efficient model often within free tier limits
# Target classification types we want to keep
TARGET_TYPES = {
    "domain", # Covers 'Functional Domain'
    "active site",
    "binding site",
    "conserved site",
    "ptm site", # Post-translational modification site
    "region",
    "repeat",
    "motif"
}

# --- TSV Output Configuration ---
TSV_OUTPUT_PATH = OUTPUT_FILE_PATH.replace('.json', '_progress.tsv')
TSV_HEADERS = ["IPR_ID", "Classification_Type", "Short_Name", "Name", "GO_Terms_Count",
               "Molecular_Function_GO_Terms", "GO_Molecular_Function_Descriptions", "Description", "Functionality_Assessment", "Explanation", "Timestamp"]

# Initialize TSV file if it doesn't exist
def initialize_tsv_file():
    if not os.path.exists(TSV_OUTPUT_PATH):
        logging.info(f"Creating new progress tracking TSV file: {TSV_OUTPUT_PATH}")
        with open(TSV_OUTPUT_PATH, 'w', newline='', encoding='utf-8') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')
            writer.writerow(TSV_HEADERS)
    else:
        logging.info(f"Found existing progress tracking TSV file: {TSV_OUTPUT_PATH}")

# Get already processed IPR IDs from TSV
def get_processed_ipr_ids():
    processed_ids = set()
    if os.path.exists(TSV_OUTPUT_PATH):
        try:
            with open(TSV_OUTPUT_PATH, 'r', encoding='utf-8') as tsvfile:
                reader = csv.reader(tsvfile, delimiter='\t')
                next(reader)  # Skip headers
                for row in reader:
                    if row and len(row) > 0:
                        processed_ids.add(row[0])  # First column is IPR_ID
            logging.info(f"Found {len(processed_ids)} already processed IPR IDs")
        except Exception as e:
            logging.error(f"Error reading existing TSV file: {e}")
    return processed_ids

# Append batch results to TSV immediately
def append_to_tsv(batch_entries, batch_results, batch_explanations):
    try:
        with open(TSV_OUTPUT_PATH, 'a', newline='', encoding='utf-8') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')
            timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            for entry in batch_entries:
                ipr_id = entry[0]
                if ipr_id in batch_results:
                    classification_type = entry[1]
                    short_name = entry[2]
                    name = entry[3]
                    go_terms_count = len(entry[4]) if entry[4] else 0

                    # Extract molecular function GO terms
                    molecular_function_go_terms = []
                    molecular_function_descriptions = []

                    if entry[4]:
                        for go_term in entry[4]:
                            # Check if GO term is molecular function (either by ID pattern or explicit category)
                            if (isinstance(go_term, list) and len(go_term) >= 3 and
                                ((go_term[0].startswith("GO:0003") or
                                 ("molecular_function" in go_term[1].lower())))):
                                # Add full GO term with ID
                                molecular_function_go_terms.append(f"{go_term[0]}: {go_term[2]}")
                                # Add just the description
                                molecular_function_descriptions.append(go_term[2])

                    # Join molecular function terms with semicolons
                    molecular_function_str = "; ".join(molecular_function_go_terms) if molecular_function_go_terms else "None"
                    molecular_function_desc_str = "; ".join(molecular_function_descriptions) if molecular_function_descriptions else "None"

                    # Use full description instead of preview
                    full_description = entry[5] if entry[5] else "None"

                    functionality = batch_results.get(ipr_id, "UNKNOWN")
                    explanation = batch_explanations.get(ipr_id, "No explanation provided")

                    writer.writerow([
                        ipr_id, classification_type, short_name, name,
                        go_terms_count, molecular_function_str, molecular_function_desc_str,
                        full_description, functionality, explanation, timestamp
                    ])
        logging.info(f"Appended {len(batch_results)} results to TSV file")
    except Exception as e:
        logging.error(f"Error appending to TSV file: {e}")

# --- Configure the Gemini API ---
logging.info(f"Configuring Gemini API with model: {MODEL_NAME}")
genai.configure(api_key=API_KEY)

generation_config = {
    "temperature": 0.1, # Very low temperature for consistent classification
    "top_p": 1,
    "top_k": 1,
    "max_output_tokens": 4096, # Reduced slightly, still ample for 10 entries
}

# Stricter safety settings might sometimes block borderline content, adjust if needed
safety_settings = [
    {"category": "HARM_CATEGORY_HARASSMENT", "threshold": "BLOCK_ONLY_HIGH"},
    {"category": "HARM_CATEGORY_HATE_SPEECH", "threshold": "BLOCK_ONLY_HIGH"},
    {"category": "HARM_CATEGORY_SEXUALLY_EXPLICIT", "threshold": "BLOCK_ONLY_HIGH"},
    {"category": "HARM_CATEGORY_DANGEROUS_CONTENT", "threshold": "BLOCK_ONLY_HIGH"},
]

model = genai.GenerativeModel(
    model_name=MODEL_NAME,
    generation_config=generation_config,
    safety_settings=safety_settings
)

# --- Load and Prepare Data ---
logging.info(f"Loading JSON data from: {JSON_FILE_PATH}...")
try:
    # Using ijson for potentially lower memory usage on large files
    # If ijson is not installed: pip install ijson
    import ijson
    unique_interpro_entries = {}
    with open(JSON_FILE_PATH, 'rb') as f: # Open in binary mode for ijson
        # Assuming top level is objects (dict), values are lists of lists
        parser = ijson.items(f, '', multiple_values=True) # Use '' for root, handle dict stream
        count = 0
        for item in parser:
            count += 1
            if count % 100 == 0:
                logging.info(f" Scanned {count} source entries...")

            # Check if the yielded item is a 2-element sequence before unpacking
            if isinstance(item, dict) and len(item) > 1:
                pro_entries = [ (source_id, entries) for source_id, entries in item.items() ] # Unpack only if it's safe
            else:
                # Log a warning if the structure is not as expected and skip this item
                logging.warning(f"Skipping unexpected item structure yielded by ijson (expected 2 elements): {type(item)}, length: {len(item)} - Content preview: {str(item)[:100]}...")
                continue # Move to the next item from the parser

            if not pro_entries:
                logging.warning(f"Skipping empty or malformed entry list for source_id {source_id}")
                continue

            entries = [e for t in pro_entries for e in t[1]]
            for entry_data in entries:
                if isinstance(entry_data, list) and len(entry_data) >= 6:
                    interpro_id = entry_data[0]
                    # Simple check if it looks like an IPR ID
                    if isinstance(interpro_id, str) and interpro_id.startswith("IPR"):
                        if interpro_id not in unique_interpro_entries:
                            # Store: IPR ID, Type, Short Name, Name, GO, Desc
                            unique_interpro_entries[interpro_id] = entry_data[0:6]
                    else:
                        logging.warning(f"Skipping entry with potentially invalid InterPro ID '{interpro_id}' for source_id {source_id}")
                else:
                    logging.warning(f"Skipping malformed entry data structure for source_id {source_id}: {entry_data}")

except FileNotFoundError:
    logging.error(f"Error: File not found at {JSON_FILE_PATH}")
    exit(1)
except ImportError:
    logging.error("Module 'ijson' not found. Please install it: pip install ijson")
    exit(1)
except Exception as e:
    logging.error(f"An error occurred during JSON loading/parsing: {e}")
    exit(1)

logging.info("JSON data loaded and parsed.")

entries_list = list(unique_interpro_entries.values())
total_entries = len(entries_list)
logging.info(f"Found {total_entries} unique InterPro entries to classify.")

# Initialize TSV file
initialize_tsv_file()

# Get already processed IPR IDs to avoid redundant processing
processed_ipr_ids = get_processed_ipr_ids()
logging.info(f"Will skip {len(processed_ipr_ids)} already processed IPR IDs")

# --- Process in Batches ---
logging.info(f"Starting processing in batches of {BATCH_SIZE}...")
final_curated_list = {} # Store all parsed IPR IDs with their binary functionality assessment
explanations_dict = {} # Store explanations for why entries are functional

# Regex to strictly parse the expected output format (now for yes/no responses)
output_pattern = re.compile(r"^(IPR\d+):\s*(yes|no)$", re.IGNORECASE)
# Regex to detect the start of explanations section
explanation_section_pattern = re.compile(r"^===\s*EXPLANATIONS\s*===\s*$", re.IGNORECASE)
# Regex to parse individual explanations
explanation_pattern = re.compile(r"^(IPR\d+):\s*(.+)$", re.IGNORECASE)

for i in range(0, total_entries, BATCH_SIZE):
    # Filter out entries that have already been processed
    batch = [entry for entry in entries_list[i:i+BATCH_SIZE] if entry[0] not in processed_ipr_ids]

    if not batch:
        logging.info(f"Skipping batch {(i // BATCH_SIZE) + 1} - all entries already processed")
        continue

    batch_number = (i // BATCH_SIZE) + 1
    start_entry_num = i + 1
    end_entry_num = min(i + BATCH_SIZE, total_entries)
    logging.info(f"--- Processing Batch {batch_number} (Entries {start_entry_num}-{end_entry_num}/{total_entries}) ---")
    logging.info(f"Batch contains {len(batch)} unprocessed entries")

    # --- Construct the Prompt ---
    prompt_parts = [
        "You are an expert bioinformatician and protein structure expert analyzing InterPro entries.",
        "Your task is to determine whether each provided InterPro entry represents a functional domain/motif/site/region of a protein.",
        "Be careful to those families of functional domains, for instance: 'IPR050224', 'Family', 'TALE_homeobox', 'Three Amino acid Loop Extension (TALE) homeobox', [], 'The TALE homeobox family is characterized by a conserved DNA-binding homeodomain' ",
        "Despite that 'IPR050224' represents a family classification, it represents a family of homoeodomain, which still counts as a functional domain/motif/site/region on proteins.",
        "Another case: ('IPR036910', 'Homologous_superfamily', 'HMG_box_dom_sf', 'High mobility group box domain superfamily', [], 'High mobility group (HMG) box domains are involved in binding DNA'), which also represents a superfamily of functional domains, which should be counted as functional domain",
        "But you have to make sure that the selected entres do not present entire protein, it must be smaller than a whole protein like a domain/motif/site, but simutaenously functional. So it can be a family of functional domains but cant be family of functional proteins",
		"Please be more stringent on determining whether a sub-region is functional, if an entry is a domain, it is not necessarily functional. You need to gave a final decision based on all the information provided by each entry including the name, the type, the description, the matched GO terms. Also remember to give concise but informative explanations on your decisions.",
        "Format your response in TWO distinct sections:",

        "SECTION 1: Functionality Assessment",
        "- One line per entry in the exact format: 'IPRXXXXXX: YES' or 'IPRXXXXXX: NO'",
        "- Answer YES if the entry represents a functional domain/motif/site/region on proteins",
        "- Answer NO if the entry does not represent a functional element (e.g., if it's merely a family classification of a protein, not a functional domain/motif/site/region)",

        "SECTION 2: Explanations",
        "- Start this section with exactly '=== EXPLANATIONS ==='",
        "- For each entry, provide one concise statement explaining your reasoning for the yes/no decision",
        "- Format each explanation as 'IPRXXXXXX: Your detailed explanation...'",

        "IMPORTANT: Keep these sections completely separate for parsing purposes.",
        "\n--- Entries for Assessment ---"
    ]

    current_batch_ids = [] # Keep track of IDs sent in this batch
    for entry in batch:
        ipr_id, classification_type, short_name, name, go_terms, description = entry
        current_batch_ids.append(ipr_id)
        go_str_list = [f"({go[0]}, {go[1]}, {go[2]})" for go in go_terms] if go_terms else ["None"]
        go_str = "; ".join(go_str_list)

        prompt_parts.append(f"\nEntry {ipr_id}:")
        prompt_parts.append(f"  InterPro ID: {ipr_id}")
        prompt_parts.append(f"  Classification Type: {classification_type}")
        prompt_parts.append(f"  Short Name: {short_name}")
        prompt_parts.append(f"  Name: {name}")
        prompt_parts.append(f"  GO Terms: {go_str}")
        prompt_parts.append(f"  Description: {description[:400]}...") # Limit description length

    prompt_parts.append("\n--- End of Entries ---")
    prompt_parts.append("\nFunctionality Assessment (Format: IPRXXXXXX: YES or IPRXXXXXX: NO):")

    final_prompt = "\n".join(prompt_parts)

    # --- Make the API Call ---
    retries = 3
    for attempt in range(retries):
        try:
            logging.info(f"Sending request to Gemini API for batch {batch_number} (Attempt {attempt + 1}/{retries})...")
            response = model.generate_content(final_prompt)

            # --- Log the full API response ---
            logging.info(f"===== FULL API RESPONSE FOR BATCH {batch_number} =====")
            if response.parts:
                response_text = response.text
                logging.info(response_text)
            else:
                logging.info("No response text received.")
            logging.info(f"===== END OF API RESPONSE FOR BATCH {batch_number} =====")

            # --- Parse and Filter the Response ---
            logging.info(f"Received response for batch {batch_number}.")
            batch_results = {}
            batch_explanations = {}
            processed_ids_in_batch = set()

            if response.parts:
                response_text = response.text
                lines = response_text.strip().split('\n')

                # Parse in two phases: binary judgments first, then explanations
                in_explanation_section = False

                for line in lines:
                    line = line.strip()
                    if not line:
                        continue # Skip empty lines

                    # Check if we're entering the explanations section
                    if explanation_section_pattern.match(line):
                        in_explanation_section = True
                        logging.info("  Found explanations section marker")
                        continue

                    if not in_explanation_section:
                        # We're in the binary judgment section
                        match = output_pattern.match(line)
                        if match:
                            parsed_ipr_id = match.group(1).upper() # Ensure IPR ID is uppercase
                            parsed_judgment = match.group(2).strip().lower() # 'yes' or 'no'
                            processed_ids_in_batch.add(parsed_ipr_id)

                            # Store the binary judgment
                            batch_results[parsed_ipr_id] = parsed_judgment
                            logging.info(f"  Parsed judgment: {parsed_ipr_id} -> {parsed_judgment}")
                        else:
                            logging.warning(f"  Could not parse judgment line or unexpected format: '{line}'")
                    else:
                        # We're in the explanations section
                        match = explanation_pattern.match(line)
                        if match:
                            parsed_ipr_id = match.group(1).upper()
                            explanation = match.group(2).strip()
                            batch_explanations[parsed_ipr_id] = explanation
                            logging.info(f"  Parsed explanation for {parsed_ipr_id} is {explanation}")
                        else:
                            # This might be a continuation of the previous explanation or an unexpected format
                            logging.warning(f"  Could not parse explanation line or unexpected format: '{line}'")

                    # Add successfully parsed results to the main collections
                    final_curated_list.update(batch_results)
                    explanations_dict.update(batch_explanations)

            else:
                # Check for safety blocks or empty response
                block_reason = response.prompt_feedback.block_reason if response.prompt_feedback else "Unknown"
                logging.warning(f"Received no response parts for batch {batch_number}. Potential safety block: {block_reason}")
                # Log safety ratings if available and needed for debugging
                if response.prompt_feedback and response.prompt_feedback.safety_ratings:
                    logging.warning("Safety ratings:")
                    for rating in response.prompt_feedback.safety_ratings:
                        logging.warning(f"   Safety Category: {rating.category}, Probability: {rating.probability}")

            # Check if all sent IDs were processed
            missing_ids = [ipr_id for ipr_id in current_batch_ids if ipr_id not in processed_ids_in_batch]
            if missing_ids:
                logging.warning(f"  Batch {batch_number}: Some IDs sent were not found in the parsed response: {missing_ids}")

            # --- After successful parsing, immediately save results to TSV ---
            if batch_results:
                append_to_tsv(batch, batch_results, batch_explanations)

            break # Break retry loop on success

        except Exception as e:
            logging.error(f"An error occurred during API call for batch {batch_number} (Attempt {attempt + 1}): {e}")
            if attempt == retries - 1:
                 logging.error(f" Batch {batch_number} failed after {retries} attempts.")
            else:
                time.sleep(5 * (attempt + 1)) # Exponential backoff before retry

    # --- Rate Limiting ---
    logging.debug("Waiting 4 second before next batch...")
    time.sleep(4.1) # Wait slightly over 1 second to be safe with rate limits

logging.info("\n--- Processing Complete ---")

# --- Save Results ---
logging.info(f"Saved {len(final_curated_list)} parsed entries with their functionality assessments.")
logging.info(f"Collected explanations for {len(explanations_dict)} entries.")

# Save with more descriptive filename
output_file_path = OUTPUT_FILE_PATH.replace('.json', '_functionality.json')
logging.info(f"Saving functionality assessments to: {output_file_path}")
try:
    with open(output_file_path, 'w') as f:
        json.dump(final_curated_list, f, indent=4)
    logging.info("Functionality assessments saved successfully.")
except Exception as e:
    logging.error(f"Failed to save functionality assessments to {output_file_path}: {e}")

# Save explanations to a separate file
explanations_output_path = OUTPUT_FILE_PATH.replace('.json', '_explanations.json')
logging.info(f"Saving explanations to: {explanations_output_path}")
try:
    with open(explanations_output_path, 'w') as f:
        json.dump(explanations_dict, f, indent=4)
    logging.info("Explanations saved successfully.")
except Exception as e:
    logging.error(f"Failed to save explanations to {explanations_output_path}: {e}")

# At the end, provide a summary of progress
total_processed = len(processed_ipr_ids) + len(final_curated_list)
total_remaining = total_entries - total_processed
logging.info(f"Progress Summary:")
logging.info(f"  Total entries: {total_entries}")
logging.info(f"  Processed entries: {total_processed} ({(total_processed/total_entries)*100:.1f}%)")
logging.info(f"  Remaining entries: {total_remaining} ({(total_remaining/total_entries)*100:.1f}%)")
logging.info(f"  Progress tracking file: {TSV_OUTPUT_PATH}")

# You can also print the final list if desired
# print("\nFinal Curated List (IPR ID: Classification):")
# for ipr_id, classification in final_curated_list.items():
#     print(f"{ipr_id}: {classification}")
