#!/usr/bin/env python3
from typing import List, Dict, Any
import pandas as pd
import plotly.graph_objects as go 
import argparse
from fasta_io import read_fasta
import sys
import os 
from orf import detect_all_frames
from frameshift import FrameshiftDetector
from report import OrfReport
from visualize import ORFS, sars_cov2_genes, display_cross_sequence_comparison
from html_report2 import generate_html_report

# Emory Foerster
# Runo Siakpebru

"""
main.py

Purpose:
    Entry point for the ORF Detection with Frameshift Identification project.

What this program will do (high level):
    1. Read nucleotide sequences from an input FASTA file.
    2. For each sequence, find all ORFs (start-to-stop codon regions).
    3. Identify the longest ORF for each sequence.
    4. Analyze whether the longest ORF suggests a potential frameshift event.
    5. Compute a simple "dominance" metric for the longest ORF.
    6. Write a summary report with results for all sequences.

Inputs:
    - A FASTA file containing one or more nucleotide sequences.

Outputs:
    - A CSV file summarizing the longest ORF and frameshift status per sequence.
"""

def create_parser() -> argparse.ArgumentParser:
    """
    Create and configure the argument parser for the simulator.

    Returns:
        A configured ArgumentParser instance with all supported
        command-line arguments.
    """
    parser = argparse.ArgumentParser(description="ORF detection with frameshift analysis.")

    parser.add_argument("-r", "--reference_fasta", type=str, default="../datasets/sequence_ref.fasta", help="Path to input Fasta file.")
    parser.add_argument("-f", "--fasta_file", type=str, help="Path to input Fasta file.")
    parser.add_argument("-m", "--min_length", type=int, default=150, help="Minimum ORF length, default = 150.")
    parser.add_argument("-oc", "--output_csv", type=str, help="File name for new output CSV file.")
    parser.add_argument("-oj", "--output_json", type=str, help="File name for new output JSON file.")
    parser.add_argument("-oh", "--output_html", type=str, help="File name for new output html file.")
    parser.add_argument("-d", "--output_dir", type=str, default="../examples", help="Directory to save output files.")
    parser.add_argument("-v", "--visualize", action="store_true", help="Shows visual ORF output in terminal.")
    return parser


def main():
    """
    Purpose:
    Will execute ORF Detection pipeline

    Inputs:
    - FASTA file that contains all sequences and parameters

    Outputs:
    - A file that contains the longest ORF sequence and frameshift status

    High-Level Steps:
    1. Parse and validate command line arguments
    2. Load all FASTA records from the input file.
    3. For each sequence, find ORFs across frames and filter by minimum length.
    4. Find the longest ORF per sequence
    5. Run frameshift analysis on the longest ORF.
    6. Collect per-sequence results into a list of dictionaries.
    7. Write results to CSV or JSON or HTML.
    """
    parser = create_parser()
    args = parser.parse_args()

    # Check FASTA argument exists
    if not args.fasta_file:
        sys.stderr.write("Error: missing FASTA file. Please enter the correct command.\n")
        return
    
    
    # Check file exists
    if not os.path.isfile(args.fasta_file):
        sys.stderr.write(f"Error: FASTA file not found for the provided path: {args.fasta_file}.\n")
        return
    
    # Check that min length is positive
    if args.min_length < 0:
         sys.stderr.write("Error Minimum length must be a positive integer.\n")

    # default to CSV if no output format provided
    if not any([args.output_csv, args.output_json, args.output_html]):
        args.output_csv = "results.csv"
        sys.stdout.write("No output format specified. Defaulting to CSV.\n")

    # Loads reference sequence and recrods for analysis
    ref_records = read_fasta(args.reference_fasta)
    records = read_fasta(args.fasta_file)
    
    # Ensures reference file only contains 1 sequence
    if len(ref_records) != 1: 
        sys.stderr.write(f"Error: reference FASTA must contain exactly 1 sequence, found {len(ref_records)}.\n")

    sys.stdout.write(f"Reference:   {ref_records[0]['ID']} ({args.reference_fasta}\n)")
    sys.stdout.write(f"Comparison:  {len(records)} sequence(s) from ({args.fasta_file})\n")
    sys.stdout.write(f"Loaded {len(records) + 1} sequences from {args.fasta_file} and {args.reference_fasta}.\n")

    records = ref_records + records

    final_results = []
    # Visualization object with SARS-CoV-2 references genes
    orfs_viz_sars = ORFS(sars_cov2_genes)
    # ORFs per sequence
    all_orfs_list = []
    # Summary metrics per sequence for comparison table
    comparison_data = []
    # Full data per sequence for HTML report
    html_records = []

    for record in records:
         seq_id = record["ID"]
         full_seq = record["Sequence"]
         seq_len = len(full_seq)

         #Detects all ORFs across all 3 reading frames above the minimum length threshold
         all_orfs = detect_all_frames(full_seq, min_length = args.min_length )
         all_orfs_list.append(all_orfs)
         # Prints colorized ORF table and gene coverage to terminal or to output file
         if args.visualize: 
            orfs_viz_sars.visualize_orf(all_orfs, seq_id, full_seq, record.get("Description", ""))
            orfs_viz_sars.display_gene_coverage(all_orfs)
         if not all_orfs:
             continue
         # Runs frameshift analysis on longest detected ORF
         detector = FrameshiftDetector(all_orfs, window=200)
         long_orf = detector.longest_orf()
         result   = detector.analyze(long_orf)
         # Prints frameshift details to terminal if frameshift detected
         if args.visualize and result["frameshift_boolean"]:
            orfs_viz_sars.display_frameshift(result, full_seq, seq_id, record.get("Description", ""))
         
         result["sequence_id"] = seq_id
         result["full_seq_length"] = seq_len   
         final_results.append(result)

         # Counts ORFS overlapping the Spike gene (diagnostic indicator for SARS-CoV-2)
         spike_orfs = sum(1 for o in all_orfs if orfs_viz_sars.get_gene_name(o["start"], o["end"], seq_id, record.get("Description", "")) == "Spike (S)")
         
         # Extract frameshift position and type for comparison table
         frameshift_details = result.get("frameshift_details", [])
         if frameshift_details:
            fs_pos = frameshift_details[0].get("shift_position", "None")
            fs_type = frameshift_details[0].get("shift_type", "None")
         else: 
            fs_pos = None
            fs_type = None
         # Appends per-sequence summary for the comparison table
         comparison_data.append({
            "sequence_id":     seq_id, 
            "total_orfs":      len(all_orfs),
            "frameshift_pos":  fs_pos,
            "frameshift_type": fs_type,
            "spike_orfs":      spike_orfs,
            })
         # Appends full record data for HTML
         html_records.append({
            "seq_id":         seq_id, 
            "all_orfs":       all_orfs,
            "result":         result,
            "full_seq":       full_seq,
            "orfs_obj":       orfs_viz_sars,
            "description":    record.get("Description", ""),
            })
         sys.stdout.write(f"Processed seq_id: {seq_id}\n")
    
    # Prints sequence comparison table to terminal 
    if args.visualize and len(comparison_data) > 1:
        display_cross_sequence_comparison(comparison_data)

    #Generates the HTML report with all sequences
    generate_html_report(html_records, output_path="../examples/report2.html")

    if final_results:
        # Pass the list of results to the engine
        report_engine = OrfReport(final_results, args.output_dir)
        
        if args.output_csv:
            report_engine.produce_report("csv")
        if args.output_json:
            report_engine.produce_report("json")
        if args.output_html:
            report_engine.produce_report("html")
    else:
        sys.stdout.write("No valid ORFs found. No reports generated.\n")


if __name__ == '__main__':
	main()
    