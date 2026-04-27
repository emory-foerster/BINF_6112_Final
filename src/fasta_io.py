#!/usr/bin/env python3

#Akshitha Ganta 

"""
fasta_io.py

Reads nucleotide sequences from a FASTA file and returns them as a list of records.
Used by main.py as the first step of the pipeline.

Input: FASTA file path
Output: list of records.

Modified by Emory: added Description field parsing from FASTA header
to support sequence annotation in visualize.py and html_report2.py.
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
                    records.append({"ID": sequence_id, "Sequence": "".join(sequence).upper(), "Description": description})
                # Removes ">" and splits header into 2 parts seq_id and description
                parts = line[1:].split(None,1) 
                sequence_id = parts[0]
                description = parts[1] if len(parts) > 1 else ""
                sequence = []
            else:
                sequence.append(line)
                
        if sequence_id is not None:
            records.append({"ID": sequence_id, "Sequence": "".join(sequence).upper(), "Description": description})

    return records 
