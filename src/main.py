#!/usr/bin/env python3

import argparse
from fasta_io import read_fasta
from orf import detect_all_frames

# Emory Foerster

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
    6. Write a summary report (CSV) with results for all sequences.

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

    parser.add_argument("-f", "--fasta_file", type=str, help="Path to input Fasta file.")
    parser.add_argument("-m", "--min_length", type=int, default=150, help="Minimum ORF length, default = 150.")
    parser.add_argument("-o", "--output_csv", type=str, help="File name for new output CSV file.")

    return parser


def main():
    """
    Purpose:
    Will execute ORF Detection pipeline

    Inputs:
    - FASTA file that contains all sequences and parameters

    Outputs:
    - A CSV file that contains the longest ORF sequence and frameshift status

    High-Level Steps:
    1. Open the FASTA file
    2. Read in all sequences
    3. Identify ORFs
    4. Find the longest ORF per sequence
    5. Analyze the potential frameshift for the longest ORF
    6. Calculate dominance metrics
    7. Write results to a CSV file
    """
    parser = create_parser()
    args = parser.parse_args()

    records = read_fasta(args.fasta_file)
    seq = records[0]["Sequence"]

    all_orfs = detect_all_frames(seq)
    meaningful_orfs = detect_all_frames(seq, min_length=args.min_length)
    print(meaningful_orfs)


if __name__ == '__main__':
	main()
    
