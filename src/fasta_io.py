#!/usr/bin/env python3

#Akshitha Ganta 

"""
fasta_io.py

Reads nucleotide sequences from a FASTA file and returns them as a list of records.
Used by main.py as the first step of the pipeline.

Input: FASTA file path
Output: list of records.
"""
def read_fasta(path):
    """
    1. Open the file
    2. For each header line (starts with '>'), start a new record.
    3. Append following sequence lines until the next header.
    4. Return all records with sequences uppercased. 
    """

    records = []
    sequence_id = None
    sequence = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
        
            if line.startswith(">"):
                if sequence_id is not None:
                    records.append({"ID": sequence_id, "Sequence": "".join(sequence).upper()})
                sequence_id = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)
                
        if sequence_id is not None:
            records.append({"ID": sequence_id, "Sequence": "".join(sequence).upper()})

    return records 
