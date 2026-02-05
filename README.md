# BINF_6112_Final
## Members
Emory Foerster | efoerste@charlotte.edu | 801504115 \
Mekhi Lucas | mlucas37@charlotte.edu | 801253831 \
Akshitha Ganta | aganta@charlotte.edu | 801290143 \
Runo Siakpebru | osiakpeb@charlotte.edu | 801233192 \
URL : https://github.com/emory-foerster/BINF_6112_Final
## License 
This project is licensed under the [GNU GPL V3 License](LICENSE) - see the [LICENSE](LICENSE) file for details.
We chose this license to make our software accessible to all, while also keeping our source code and any other changes under the same license. 


ORF Detection with Frameshift Identification and Analysis 

Extend ORF analysis by incorporating frameshift detection.

Idea:
Build a python tool that scans nucleotide FASTA sequences to identify ORF frames and detect candidates
that show a likely frameshift detection. 

Core functionality:
- Read one or more nucleotide sequences from FASTA.
- Find all ORFs across reading frames.
- Identify the longest ORF in each sequence
- Detect whether the ORF includes a potential **frameshift event**
  
Report:
- The longest ORF containing a frameshift, if present, or
- The longest ORF without a frameshift, if no frameshift is detected.
- Dominance of the longest ORF compared to total ORF in sequence  

Outcome:
- A JSON summary per sequence showing the results.
- The project should explain the logic used to detect frameshifts and discuss limitations of the approach
- Find the coverage of the longest ORF with a frameshift using bam
- Using sequences in a family do a similarity score between the longest ORF's with a frameshift
- Compare longest ORF to a  database to find closely relative proteins

Project Structure:
BINF_6112_FINAL/
  README.md
  main.py
  fasta_io.py
  orf.py
  frameshift.py
  report.py

