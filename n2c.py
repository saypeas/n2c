## This version is identical to Fabricio Kury's script but is commented by Claude.

"""
NDC to ATC Classification Mapping Script

OVERVIEW:
This script processes a list of National Drug Codes (NDCs) and maps them to 
Anatomical Therapeutic Chemical (ATC) classification codes using the RxNorm 
and RxClass APIs from the National Library of Medicine.

FUNCTIONALITY:
- Reads NDC codes from a text file (one per line)
- Maps NDCs to ATC-4 or ATC-5 classification codes through RxNorm CUIs
- Supports two mapping types:
  * ATC-4: Maps directly from RxCUI to ATC-4 classes
  * ATC-5: Maps from RxCUI → ingredients → ATC-5 classes
- Uses persistent caching to avoid redundant API calls
- Provides progress tracking with time estimates
- Handles interruption gracefully (press 'P' to stop)
- Outputs results to CSV format

WORKFLOW:
1. NDC → RxCUI (RxNorm Concept Unique Identifier)
2a. For ATC-4: RxCUI → ATC-4 classes
2b. For ATC-5: RxCUI → ingredient CUIs → ATC-5 classes
3. Save all mappings to CSV file

DEPENDENCIES:
- requests: HTTP API calls
- csv: CSV file handling
- shelve: Persistent caching
- keyboard: Interrupt detection
- argparse: Command line argument parsing

USAGE:
python script.py input_file.txt --mapping atc4 [--output_file output.csv]
python script.py input_file.txt --mapping atc5 [--output_file output.csv]
"""

import argparse
import requests
import csv
import time
import os
import shelve
import keyboard
from collections import defaultdict
from datetime import datetime

# Constants for API URLs
RXNORM_API_URL = "https://rxnav.nlm.nih.gov/REST/rxcui"
RXCLASS_API_URL = "https://rxnav.nlm.nih.gov/REST/rxclass/class/byRxcui"

def validate_ndc(ndc):
    """
    Validate the format of the NDC (National Drug Code).
    
    NDC codes should be 10-12 digits long and may contain hyphens.
    Valid formats: XXXXXXXXXX, XXXXX-XXX-XX, XXXX-XXXX-XX, etc.
    
    Args:
        ndc (str): The NDC code to validate
        
    Returns:
        bool: True if valid format, False otherwise
    """
    if len(ndc) in [10, 11, 12] and ndc.replace("-", "").isdigit():
        return True
    print(f"Invalid NDC format: {ndc}")
    return False

def format_time(seconds):
    """
    Format time in seconds to a human-readable string.
    
    Converts seconds into hours, minutes, and seconds format for 
    better progress reporting.
    
    Args:
        seconds (int): Time in seconds
        
    Returns:
        str: Formatted time string (e.g., "2 hours, 30 minutes, 15 seconds")
    """
    periods = [
        ('hour', 3600),
        ('minute', 60),
        ('second', 1)
    ]
    time_str = []
    for period_name, period_seconds in periods:
        if seconds >= period_seconds:
            period_value, seconds = divmod(seconds, period_seconds)
            if period_value > 1:
                period_name += 's'
            time_str.append(f"{period_value} {period_name}")
    return ', '.join(time_str)

def get_rxcui_from_ndc(ndc, cache):
    """
    Convert NDC (National Drug Code) to RxCUI (RxNorm Concept Unique Identifier).
    
    Makes API call to RxNorm to get the unique identifier for the drug.
    Uses caching to avoid redundant API calls for the same NDC.
    
    Args:
        ndc (str): The NDC code to look up
        cache (shelve.Shelf): Persistent cache for storing results
        
    Returns:
        str or None: RxCUI if found, None if not found or error occurred
    """
    cache_key = f"ndc:{ndc}"
    if cache_key in cache:
        return cache[cache_key]
    try:
        response = requests.get(f"{RXNORM_API_URL}?idtype=NDC&id={ndc}", headers={"Accept": "application/json"})
        response.raise_for_status()
        data = response.json()
        rxcui = data.get('idGroup', {}).get('rxnormId', [None])[0]
        cache[cache_key] = rxcui
        return rxcui
    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
    except requests.exceptions.RequestException as req_err:
        print(f"Request error occurred: {req_err}")
    except ValueError as json_err:
        print(f"JSON decode error: {json_err}")
    cache[cache_key] = None
    return None

def get_atc4_classes_from_rxcui(rxcui, cache):
    """
    Get ATC-4 classification codes from RxCUI using RxClass API.
    
    ATC-4 codes represent therapeutic subgroups (4th level of ATC hierarchy).
    Makes API call to RxClass to get all ATC classifications for the drug.
    
    Args:
        rxcui (str): RxNorm Concept Unique Identifier
        cache (shelve.Shelf): Persistent cache for storing results
        
    Returns:
        list: List of ATC-4 class codes, empty list if none found or error
    """
    cache_key = f"atc4:{rxcui}"
    if cache_key in cache:
        return cache[cache_key]
    try:
        response = requests.get(f"{RXCLASS_API_URL}?rxcui={rxcui}&classTypes=ATC", headers={"Accept": "application/json"})
        response.raise_for_status()
        data = response.json()
        classes = data.get('rxclassDrugInfoList', {}).get('rxclassDrugInfo', [])
        atc_classes = [cls['rxclassMinConceptItem']['classId'] for cls in classes if cls['rxclassMinConceptItem']['classType'] == 'ATC1-4']
        cache[cache_key] = atc_classes
        return atc_classes
    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
    except requests.exceptions.RequestException as req_err:
        print(f"Request error occurred: {req_err}")
    except ValueError as json_err:
        print(f"JSON decode error: {json_err}")
    cache[cache_key] = []
    return []

def get_ingredients_from_rxcui(rxcui, cache):
    """
    Get ingredient RxCUIs from a drug's RxCUI.
    
    Retrieves all active ingredients (TTY = 'IN') for a given drug.
    This is needed for ATC-5 mapping since ATC-5 codes are typically
    associated with individual ingredients rather than drug products.
    
    Args:
        rxcui (str): RxNorm Concept Unique Identifier for the drug
        cache (shelve.Shelf): Persistent cache for storing results
        
    Returns:
        list: List of ingredient RxCUIs, empty list if none found or error
    """
    cache_key = f"ingredients:{rxcui}"
    if cache_key in cache:
        return cache[cache_key]
    try:
        url = f"https://rxnav.nlm.nih.gov/REST/rxcui/{rxcui}/allrelated.json"
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        ingredients = []
        for concept in data['allRelatedGroup']['conceptGroup']:
            if concept.get('tty') == 'IN':
                ingredients.extend([ing['rxcui'] for ing in concept['conceptProperties']])
        cache[cache_key] = ingredients
        return ingredients
    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
    except requests.exceptions.RequestException as req_err:
        print(f"Request error occurred: {req_err}")
    except ValueError as json_err:
        print(f"JSON decode error: {json_err}")
    cache[cache_key] = []
    return []

def get_atc5_classes_from_ingredient(ingredient_cui, cache):
    """
    Get ATC-5 classification codes from an ingredient RxCUI.
    
    ATC-5 codes represent specific chemical substances (5th level of ATC hierarchy).
    These are typically mapped to individual ingredients rather than drug products.
    
    Args:
        ingredient_cui (str): RxCUI for the drug ingredient
        cache (shelve.Shelf): Persistent cache for storing results
        
    Returns:
        list: List of ATC-5 class codes, empty list if none found or error
    """
    cache_key = f"atc5:{ingredient_cui}"
    if cache_key in cache:
        return cache[cache_key]
    try:
        url = f'https://rxnav.nlm.nih.gov/REST/rxcui/{ingredient_cui}/property.json?propName=ATC'
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        atc_classes = []
        # Check if 'propConceptGroup' and 'propConcept' exist in the response
        if 'propConceptGroup' in data and 'propConcept' in data['propConceptGroup']:
            for prop in data['propConceptGroup']['propConcept']:
                if prop['propName'].startswith('ATC'):
                    atc_classes.append(prop['propValue'])
        cache[cache_key] = atc_classes
        return atc_classes
    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
    except requests.exceptions.RequestException as req_err:
        print(f"Request error occurred: {req_err}")
    except ValueError as json_err:
        print(f"JSON decode error: {json_err}")
    cache[cache_key] = []
    return []

def generate_output_filename(input_file, mapping_type):
    """
    Generate output filename based on input filename and mapping type.
    
    Creates a descriptive filename that indicates the type of mapping performed.
    
    Args:
        input_file (str): Path to the input file
        mapping_type (str): Type of mapping ('atc4' or 'atc5')
        
    Returns:
        str: Generated output filename (e.g., "input_ATC4_classes.csv")
    """
    base, ext = os.path.splitext(input_file)
    return f"{base}_{mapping_type.upper()}_classes.csv"

def generate_cache_filename(input_file, mapping_type):
    """
    Generate cache filename for persistent storage of API results.
    
    Creates a cache file specific to the input file and mapping type
    to avoid cache conflicts between different runs.
    
    Args:
        input_file (str): Path to the input file
        mapping_type (str): Type of mapping ('atc4' or 'atc5')
        
    Returns:
        str: Generated cache filename (e.g., "input_atc4_cache.shelve")
    """
    base = os.path.splitext(input_file)[0]
    return f"{base}_{mapping_type}_cache.shelve"

def process_ndc_list(input_file, output_file, mapping_type, cache):
    """
    Main processing function that maps NDCs to ATC codes.
    
    Reads NDCs from input file, performs the specified mapping (ATC-4 or ATC-5),
    and writes results to CSV. Handles progress tracking, caching, and graceful
    interruption.
    
    Process flow:
    1. Read and validate NDCs from input file
    2. For each NDC:
       - Get RxCUI from NDC
       - If ATC-4: Get ATC-4 classes directly from RxCUI
       - If ATC-5: Get ingredients from RxCUI, then ATC-5 classes from ingredients
    3. Save results to CSV with deduplication
    
    Args:
        input_file (str): Path to file containing NDCs (one per line)
        output_file (str): Path for output CSV file
        mapping_type (str): Either 'atc4' or 'atc5'
        cache (shelve.Shelf): Persistent cache for API results
    """
    ndcs_with_atc = 0
    results = []
    ndcs_no_rxcui = 0

    # Read NDCs from input file
    with open(input_file, 'r') as file:
        ndcs = [line.strip() for line in file if validate_ndc(line.strip())]

    total_ndcs = len(ndcs)
    api_calls = 0
    start_time = datetime.now()

    try:
        for idx, ndc in enumerate(ndcs, start=1):
            # Check for keyboard interrupt (P key)
            if keyboard.is_pressed('p'):
                raise KeyboardInterrupt

            # Get RxCUI from NDC
            rxcui = get_rxcui_from_ndc(ndc, cache)

            if rxcui:
                if mapping_type == 'atc4':
                    # Get ATC-4 classes from RxCUI
                    atc_classes = get_atc4_classes_from_rxcui(rxcui, cache)
                    if atc_classes:
                        ndcs_with_atc += 1
                    for atc_class in atc_classes:
                        results.append({'NDC': ndc, 'ATC4 Class': atc_class})
                elif mapping_type == 'atc5':
                    # Get ingredients from RxCUI
                    ingredients = get_ingredients_from_rxcui(rxcui, cache)
                    atc5_classes = set()
                    for ing_cui in ingredients:
                        # Get ATC-5 classes from ingredient
                        classes = get_atc5_classes_from_ingredient(ing_cui, cache)
                        atc5_classes.update(classes)
                    if atc5_classes:
                        ndcs_with_atc += 1
                        for atc_class in atc5_classes:
                            results.append({'NDC': ndc, 'ATC5 Class': atc_class})
                    else:
                        results.append({'NDC': ndc, 'ATC5 Class': 'No ATC Mapping Found'})
                else:
                    print(f"Unknown mapping type: {mapping_type}")
                    return
            else:
                ndcs_no_rxcui += 1
                if mapping_type == 'atc5':
                    results.append({'NDC': ndc, 'ATC5 Class': 'No RxCUI Found'})
                # For 'atc4', we simply skip NDCs without RxCUI (as per original script)

            # Save the cache every 10 iterations to minimize data loss
            if idx % 10 == 0:
                cache.sync()

            api_calls += 1
            # Update progress every 10 NDCs
            if idx % 10 == 0 or idx == total_ndcs:
                current_time = datetime.now()
                elapsed_time = current_time - start_time
                estimated_total_time = (elapsed_time / idx) * total_ndcs
                remaining_time = estimated_total_time - elapsed_time
                remaining_seconds = int(remaining_time.total_seconds())
                print(f"Processing {idx}/{total_ndcs}: {ndc} - Estimated completion in {format_time(remaining_seconds)}")

    except KeyboardInterrupt:
        completion_percentage = (ndcs_with_atc / idx) * 100
        print(f"Interrupted. {completion_percentage:.2f}% of NDCs queried have at least one ATC class associated.")
        cache.sync()

    # De-duplicate the results if necessary
    unique_results = {tuple(result.items()) for result in results}
    unique_results = [dict(result) for result in unique_results]

    # Write results to CSV
    with open(output_file, 'w', newline='') as csvfile:
        if mapping_type == 'atc4':
            fieldnames = ['NDC', 'ATC4 Class']
        elif mapping_type == 'atc5':
            fieldnames = ['NDC', 'ATC5 Class']
        else:
            print(f"Unknown mapping type: {mapping_type}")
            return
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(unique_results)

    print(f"Results written to {output_file}")
    completion_percentage = (ndcs_with_atc / total_ndcs) * 100
    print(f"Completed: {completion_percentage:.2f}% of NDCs have at least one ATC class associated.")

if __name__ == "__main__":
    """
    Command-line interface for the NDC to ATC mapping script.
    
    Parses command-line arguments and orchestrates the mapping process.
    Handles file naming and cache management.
    """
    parser = argparse.ArgumentParser(description='Process a list of NDCs and map them to ATC classes.')
    parser.add_argument('input_file', type=str, help='Input file containing NDCs (one per line).')
    parser.add_argument('--mapping', type=str, choices=['atc4', 'atc5'], required=True,
                        help='Type of mapping to perform: atc4 or atc5.')
    parser.add_argument('--output_file', type=str, help='Output file to write results to.')

    args = parser.parse_args()

    if not args.output_file:
        output_file = generate_output_filename(args.input_file, args.mapping)
    else:
        output_file = args.output_file

    cache_file = generate_cache_filename(args.input_file, args.mapping)

    print("Will start querying the NDCs. If you need to interrupt, press P.")

    # Use shelve to create a persistent cache in the same directory as the input file
    with shelve.open(cache_file) as cache:
        process_ndc_list(args.input_file, output_file, args.mapping, cache)
