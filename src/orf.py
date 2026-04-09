#!/usr/bin/env python3

#Emory Foerster 

from fasta_io import read_fasta

def parse_codons(seq: str, frame:int) -> list[tuple[str, int]]:
    """
    Purpose:
        Split a nucleotide sequence into codons starting
        given fram offset. 
    Input:
        sequence (str): nucleotide sequence list from fastaio.py.
        frame (int): reading frame offset -0 , 1 , or 2
    Output:
        A list of tuples: (codon_string, start_position)
        Example: [("ATG", 0), ("TCA", 3), ("TAA", 6)]
    High‑level steps:
        1. start index at frame 
        2. scan sequence 3 nucleotides at a time
            - only includes complete codons
        3. records each codon and its start index in the sequence
        4. returns list of tuples with codon + start index [("codon", start_index), ]
    Raise Errors: 
        ValueError: if sequence is empty or not a string type
        ValueError: if frame is not 0,1 or 2
    """
    if not seq or not isinstance(seq, str):
        raise ValueError("Sequence must be a non-empty and of string type.")
    if frame not in (0,1,2):
        raise ValueError(f"Frame must be 0, 1, or 2. Frame is {frame}.")

    codons = []
    for i in range(frame, len(seq) -2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3:
            codons.append((codon,i))
    return codons

def detect_ORF(codon_list: list[tuple[str,int]], frame: int) -> list[dict]:
    """
    Purpose:
        To find all ORFs in a sequence using 3 frames and saving them into dictionaries. 
    Input:
        sequence (str): nucleotide sequence list from fastaio.py.
    Output:
        A dictionary with the reading frames and orfs cooresponding to them
        - seq      (str): nulecotide string of the ORF (start to stop condons inclusive)
        - start    (int): start positon of ORF
        - end      (int): ending position 
        - length   (int): length of ORF
        - frame    (int): reading frame 0,1,2
    High‑level steps:
        1. Read in sequence 
        2. For each reading frames (0,1,2)
            a. initialize tracking variable to store data
            b. scan sequence through each reading frame one at a time
            c. scan sequence 3 nucleotides at a time
            d. start recording sequence when hit ATG codon 
                - if start codon comes before stop, create new ORF and record nucleotides in both
            e. stop recording sequence when hit stop codons

        3. collect all ORF from reading frames
        4. return list of ORF and the additional information 
    """
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []
    recording = False
    current_orf = ""
    start_pos = None

    for codon, position in codon_list:
        if not recording:
            if codon == "ATG":
                recording = True
                start_pos = position
                current_orf = codon
        
        else: 
            if codon in stop_codons:
                current_orf += codon
                end_pos = position + 3
                length = end_pos - start_pos
                orf = {"seq": current_orf, "start":start_pos, "end":end_pos, "length":length, "frame":frame}
                orfs.append(orf)
                recording = False
                current_orf = ""
            else: 
                current_orf += codon
    return orfs
    
def detect_all_frames(seq: str, min_length: int = 0) -> list[dict]:
    """
    Purpose:
        To run detect_ORF across all three reading frames in a sequence and returns every ORF found. 
    Input:
        sequence (str): nucleotide sequence (uppercase)
    Output:
        A list that combines all ORF dictionaries created by detect_ORF 
    High‑level steps:
        1. Run detect_ORF for frame 0,1,2
        2. Combine list of dictionaries across all three frames into one list
        3. Return the combines list
    """
    ORFs = []
    for i in range(3):
        codon_list = parse_codons(seq,i)
        orf = detect_ORF(codon_list, i)
        if min_length > 0:
            orf = [o for o in orf if o['length']>= min_length]
        ORFs.extend(orf)
    return ORFs

def main():
    record = read_fasta("../datasets/Covid_GCF_009858895.2/sample.fasta")
    seq = record[0]["Sequence"]
    all_orfs = detect_all_frames(seq)
    meaningful_orfs = detect_all_frames(seq, min_length=150)
    print(all_orfs)
    


if __name__ == '__main__':
    main()



