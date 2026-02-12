#!/usr/bin/env python3
"""
report.py

Writes the final results to an output file.
Supports CSV or JSON output, based on the user's chosen format.
Used as the last step of the pipeline called by main.py.

# RESEARCH: 
We want this script to have an output that visualizes where the frameshifts occur and give a summary of the location
and type of possible frameshift. 
- research areas: how to produced visualize and formatted output. 
"""

def produce_report(results, output_path, output_format):
    """
    Purpose:
        Write results to disk in the requested format.

    Input:
        results (list[dict]): one dict per sequence, including ORF and frameshift fields.
        output_path (str): filename to write.
        output_format (str): "json" or "csv".

    Output:
        None (writes a file).

    High-level steps:
        1. Normalize output_format to lowercase.
        2. If output_format == "json":
        - write results as a JSON list of dicts
        3. If output_format == "csv":
        - write header row using dict keys
        - write one row per result
        4. Otherwise:
        - raise ValueError for unsupported format.
    """
    pass
