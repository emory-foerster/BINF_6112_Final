#!/usr/bin/env python3
#Runo Siakpebru
from typing import List, Dict, Any
import pandas as pd
import os
import plotly.graph_objects as go 
import sys


class OrfReport:
    """
    Generate output files (JSON, CSV, HTML) from ORF + frameshift results.

    Purpose:
        Take analysis results and write them to one or more report formats.

    Input:
        results (List[Dict[str, Any]]): list of per-sequence result dictionaries.
        output_dir (str): directory where report files will be saved.

    Output:
        Report files produced
    """

    def __init__(self, results: List[Dict[str, Any]], output_dir: str = ".") -> None:
        """
        Initialize the report generator.

        Purpose: Initializes the report with raw result data.

        Inputs:
            results: list of per-sequence result dictionaries.
            output_dir: folder to write report files into.

        High-level steps:
            1. Store results in self.results.
            2. Store output_dir in self.output_dir.
            3. Create output_dir if it does not exist.
        """
        self.results = results
        self.output_dir = output_dir

        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir, exist_ok=True)

    def write_json(self, filename: str = "report.json") -> None:
        """
        Purpose: Write results to a JSON file.

        Input:
           file name(str): file name default report.json.

        Output: JSON file.

        High-level steps:
            1. If results is empty, warn and stop.
            2. append the filename to self.output_dir
            3. Write the results list to JSON using pandas. 
    """
        if not self.results:
            sys.stderr.write("Warning: No results to write to JSON.")
            return

        output_path = os.path.join(self.output_dir, filename)
        
        df= pd.DataFrame(self.results)
        df.to_json(output_path, orient='records', indent = 4)
        print(f"JSON report successfully saved to {output_path}")
        


    def write_csv(self, filename: str = "report.csv") -> None:
        """
        Purpose: Write results to a CSV file.

        Input:
           file name(str): file name default report.csv.

        Output: CSV file.

        High-level steps:
            1. If results is empty, warn and stop.
            2. append the filename to self.output_dir
            3. Extract the info from self.results and build the structure of
            4. Write the results list to JSON using pandas.     
        """
        if not self.results:
            sys.stderr.write("No results to write to CSV.\n")
            return
        
        output_path = os.path.join(self.output_dir, filename)
        
        rows = []
        for res in self.results:
            # Extract primary shift info if it exists
            details = res.get("frameshift_details", [])
            primary = details[0] if details else {}
            neighbor = primary.get("neighboring_orf", {})

            row = {
                "sequence_id": res.get("sequence_id"),
                "frameshift_suspected": res.get("frameshift_boolean"),
                "frameshift_position": primary.get("shift_position"),
                "shift_type": primary.get("shift_type"),
                "dominance": res.get("dominance_ratio"),
                "longest_orf_start": res.get("start"),
                "longest_orf_end": res.get("end"),
                "longest_orf_length_nt": res.get("length"),
                "overlap_length": neighbor.get("length")
            }
            rows.append(row)

        df= pd.DataFrame(rows)

        cols_order = [
            "sequence_id", "frameshift_suspected", "frameshift_position", 
            "shift_type", "dominance", "longest_orf_start", 
            "longest_orf_end", "longest_orf_length_nt", "overlap_length"
        ]

        cols = [c for c in cols_order if c in df.columns]
        
        df[cols].to_csv(output_path, index=False)
        print(f"CSV report successfully saved to {output_path}")


    def write_html(self, filename: str = "report.html"):
        """
        Purpose:
            Write an HTML report for all sequences. 
        Args:
            results: List of dictionaries containing ORF and frameshift data.
            output_path: The file path to save the HTML dashboard.
        """
        if not self.results:
            return

        fig = go.Figure()

        # Iterate through each processed sequence
        for res in self.results:
            seq_id = res.get("sequence_id", "Unknown")
            full_len = res.get("full_seq_length", res.get("end", 0) + 100)
            
            # 1. Draw the Full Sequence Backbone
            fig.add_trace(go.Bar(
                x=[full_len],
                y=[seq_id],
                orientation='h',
                marker_color='rgba(211, 211, 211, 0.3)', # Light Gray
                hoverinfo='skip',
                showlegend=False
            ))

            # 2. Draw the Longest ORF (The primary gene)
            fig.add_trace(go.Bar(
                x=[res.get("length")],
                y=[seq_id],
                base=res.get("start"),
                orientation='h',
                name=f"Main ORF ({seq_id})",
                marker_color='mediumseagreen',
                hovertemplate=(
                    f"<b>Primary ORF</b><br>"
                    f"Coords: {res.get('start')} - {res.get('end')}<br>"
                    f"Length: {res.get('length')}nt<extra></extra>"
                )
            ))


        output_path = os.path.join(self.output_dir, filename)
        fig.write_html(output_path)
        print(f"Interactive HTML dashboard successfully saved to {output_path}")


    def produce_report(self, output_format: str) -> None:
        """
        Purpose: Main entry point used by main.py to write outputs.

        Input: output_format (str): "json", "csv", or "html".

        Output: writes a file

        High-level steps:
            - If format is json -> call write_json().
            - If format is csv -> call write_csv().
            - If format is html -> call write_html().
            - Otherwise raise ValueError for unsupported format.
        """
        fmt = output_format.strip().lower()

        if fmt == "json":
            self.write_json()
        elif fmt == "csv":
            self.write_csv()
        elif fmt == "html":
            self.write_html()
        else:
            sys.stderr.write(f"Unsupported report format: '{output_format}'.Please use 'json', 'csv', or 'html'.")

#report1= OrfReport(results,"../examples")
#report1.produce_report('csv')