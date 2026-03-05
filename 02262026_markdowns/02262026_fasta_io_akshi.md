# fasta_io.py (Work Summary)

## What it does
Reads a FASTA file and converts it into a list of sequence records for the rest of the pipeline.

## Input
- `path` (str): path to a FASTA file containing one or more nucleotide sequences

## Output
- `records` (list[dict]): each record is a dictionary like:
  - `"ID"`: sequence ID (first token after `>`)
  - `"Sequence"`: full nucleotide sequence as one uppercase string

## Key behavior
- Reads the file line by line (works for large files)
- Combines multi-line FASTA sequences into one string per entry
- Uppercases sequences before returning

