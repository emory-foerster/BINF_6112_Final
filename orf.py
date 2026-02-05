#!/usr/bin/env python3

def detect_ORF(seq,frame):
     """
     Purpose:
         To find all ORFs in a sequence using 3 frames and saving them into dictionaries. 
     Input:
         sequence (str): nucleotide sequence list from fastaio.py.
     Output:
         A dictionary with the reading frames and orfs cooresponding to them
         - start: start positon of ORF
         - end: ending position 
         - length: length of ORF
         - frame: reading frame 0,1,2
     High‑level steps:
         1. Read in sequence 
         2. For each reading frames (0,1,2)
           a. initialize tracking variable to store data
           b. scan sequence through each reading frame one at a time
           c. scan sequence 3 nucleotides at a time
           d. start recording sequence when hit ATG codon 
               - if start codon comes before stop, create new ORF and record nucleotides in both
           e.stop recording sequence when hit stop codons

         3. collect all ORF from reading frames
         4. return list of ORF and the additional information 
     """
     pass


def main():
     detect_ORF()

if __name__ == '__main__':
	  main()


