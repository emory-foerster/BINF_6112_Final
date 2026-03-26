#!/usr/bin/env python3
#Runo Siakpebru
from typing import List, Dict, Any
import pandas as pd
import os
results = [
    {
        "sequence_id": "seqA",
        "frameshift_suspected": True,
        "frameshift_position": 600,
        "shift_type": "+1",
        "dominance": 0.78,
        "longest_orf_start": 120,
        "longest_orf_end": 1320,
        "longest_orf_length_nt": 1200
    }
]

"""
report.py

Writes the final results to an output file.
Supports CSV or JSON output, based on the user's chosen format.
Used as the last step of the pipeline called by main.py.

# RESEARCH: 
We want this script to have an output that visualizes where the frameshifts occur and give a summary of the location
and type of possible frameshift. 
- research areas: how to produced visualize and formatted output. 

Ideas from research: good article https://www.biorxiv.org/content/10.1101/2021.06.10.447953v2.full

I think we should use Plotly and Dash to create an html as an output which displays:
-type of frameshift (+1, +2) vs. length of ORF (similar to Figure 2)
-Longest ORF proportion of the whole sequence (visual and/or percentage)
-A frameshift bar with highlighted areas of frameshift potential spots(also starting their position)
-display in someway the read depth and coverage of longest ORF
-JSON summary

"""

def validate_results(results):
    """
    Purpose:
        Ensure results is in the expected shape before writing any output.
    High-level steps:
        1. Check results is a list.
        2. Check each item is a dictionary.
        3. Optionally check required keys exist (sequence_id, longest_orf fields, frameshift fields).
        4. If checks fail, raise ValueError with a clear message.
    """
    pass 

def write_json(results: List[Dict[str, Any]], output_path: str) -> None:
    """
    Write results to a JSON file.

    Args:
        results: A list of dictionaries, one per sequence.
        output_path: The file path to save the JSON.
    
    """
    if not results:
        print("Warning: No results to write to JSON.")
        return

    # If provided a directory, append a default filename
    if os.path.isdir(output_path):
        output_path = os.path.join(output_path, "report.json")
        print(f"Directory detected. Saving as: {output_path}")
    
    df= pd.DataFrame(results)
    df.to_json(output_path, orient='records', indent = 4)
    print(f"JSON report successfully saved to {output_path}")
        


def write_csv(results, output_path):
    """
     Write results to a CSV file.

    results: A list of dictionaries, one per sequence. output_path: The file path to save the CSV.      
    """
    if not results:
        print("Warning: No results to write to CSV.")
        return
    
    # If provided a directory, append a default filename
    if os.path.isdir(output_path):
        output_path = os.path.join(output_path, "report.csv")
        print(f"Directory detected. Saving as: {output_path}")

    df= pd.DataFrame(results)
    
    cols_order = [
        "sequence_id", 
        "frameshift_suspected", 
        "frameshift_position", 
        "shift_type", 
        "dominance", 
        "longest_orf_start", 
        "longest_orf_end", 
        "longest_orf_length_nt"
    ]

    cols = list()
    for c in cols_order:
        if c in df.columns:
            cols.append(c)
    
    df[cols].to_csv(output_path, index=False)
    print(f"CSV report successfully saved to {output_path}")





def write_html(results, output_path):
    """
    Purpose:
        Write a lightweight HTML report that includes a table and simple visualization bars.

    Input:
        results (list[dict]): one dict per sequence.
        output_path (str): filename to write.

    Output:
        None (writes a file)

    High-level steps:
        1. Build an HTML header with basic CSS styling.
        2. For each result dict:
            a. Extract key fields (sequence_id, lengths, frameshift info, dominance).
            b. Generate a frameshift bar using make_frameshift_bar().
            c. Add a table row to the HTML.
        3. Close the HTML document and write it to output_path.
    """
    pass


def produce_report(results, output_path, output_format):
    """
    Purpose:
        Main entry point used by main.py to write outputs.

    Input:
        output_format (str): "json", "csv", or "html".

    Output: writes a file

    High-level steps:
        - If format is json -> call write_json().
        - If format is csv -> call write_csv().
        - If format is html -> call write_html().
        - Otherwise raise ValueError for unsupported format.
    """
    fmt = output_format.strip().lower()

    if fmt == "json":
        write_json(results, output_path)
    elif fmt == "csv":
        write_csv(results, output_path)
    else:
        raise ValueError("Only 'json' and 'csv' options available now.")
produce_report(results, "report.csv", "csv")