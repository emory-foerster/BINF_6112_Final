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