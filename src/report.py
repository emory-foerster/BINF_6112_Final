#!/usr/bin/env python3
import json
import csv
#Runo Siakpebru

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
    if not isinstance(results, list):
        raise ValueError("results must be a list")

    if len(results) == 0:
        raise ValueError("results is empty")

    for i, item in enumerate(results):
        if not isinstance(item, dict):
            raise ValueError(f"results[{i}] must be a dict")


def normalize_format(output_format):
    """
    Purpose:
        Normalize the output format string so the rest of the code is consistent.

    Input:
        output_format (str): user format choice like "json", "csv", "html".

    Output:
        str: normalized format (lowercase, trimmed).

    High-level steps:
        1. Convert output_format to string, strip whitespace, lowercase it.
        2. If empty, raise ValueError.
        3. Return normalized value.
    """
    pass


def write_json(results, output_path):
    """
    Purpose:
        Write results to a JSON file.

    Input: results (list[dict]): one dict per sequence. output_path (str): filename to write.

    Output: writes a file

    High-level steps:
        1. Open output_path for writing.
        2. Write results as a JSON list (pretty print optional).
        3. Close the file.
    """
    validate_results(results)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
        


def write_csv(results, output_path):
    """
    Purpose:
        Write results to a CSV file.

    Input:
        results (list[dict]): one dict per sequence.
        output_path (str): filename to write.

    Output:
        None (writes a file)
        
    """
    fieldnames = []
    s = set()

    for row in results:
        for key in row.keys():
            if key not in s:
                s.add(key)
                fieldnames.append(key)

    # Write CSV
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()

        for row in results:
            writer.writerow(row)





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