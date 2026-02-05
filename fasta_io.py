#!/usr/bin/env python3
"""
fasta_io.py

Reads nucleotide sequences from a FASTA file and returns them as a list of records.
Used by main.py as the first step of the pipeline.

Input: FASTA file path
Output: list of records.
"""
def read_fasta():
    """
    1. Open the file
    2. For each header line (starts with '>'), start a new record.
    3. Append following sequence lines until the next header.
    4. Return all records with sequences uppercased. 
    """
    pass
