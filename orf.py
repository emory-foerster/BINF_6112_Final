#!/usr/bin/env python3

#Emory Foerster 

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
    """
    codons = []
    for i in range(frame, len(seq) -2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3:
            codons.append((codon,i))
    return codons

def detect_ORF(seq: str, frame: int) -> list[dict]:
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
    pass

def main():
    seq = "ATGCGTACGTTAGCTAGCTAGCTAGCTA"
    print(parse_codons(seq, 1))


if __name__ == '__main__':
    main()



