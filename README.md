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

## Project
ORF Detection with Frameshift Identification and Analysis 

Extend ORF analysis by incorporating frameshift detection.

## Idea:
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

## Requirements:
- Python V3
- Conda binf_final environment from dependencies.txt
   - Using conda create -n binf_final -f dependencies.txt


## Project Structure:
```
ORF Frameshift Analysis & Comparison/
├── README.md
├── LICENSE
├── ENVIRONMENT.md
├── environment.yml
├── run_test.sh
├── pseudocode/
│   ├── pseudocode.txt
│   └── flowcharts.txt
├── datasets/
│   ├── sequences.fasta
│   └── reference_sequence.fasta
├── src/
│   ├── __init__.py
│   ├── main.py
│   ├── fasta_io.py
│   ├── orf.py
│   ├── frameshift.py
│   ├── report.py
│   ├── visualize.py
│   ├── WIP_report_2.py
│   └── html_report2.py
└── out_results/
    ├── frame_job_err.txt
    ├── frame_job_out.txt
    ├── frameshift_plot.html
    ├── report.csv
    ├── report.json
    ├── results.html
    └──test_report.html
```

1. fasta_io.py
   - input: FASTA file
   - output: Dictionary of sequences
2. orf.py
   - input: fasta_io.py (Dictionary of sequences)
   - output: List of longest ORF for each entry
3. frameshift.py
   - input: orf.py (list of longest ORF's for each entry)
   - output: tbd (we are researching) - some type of text file
4. report.py
   - input: frameshift.py (tbd)
   - output: file - visualizing the potential frameshifts
5. main.py
   - input: fasta_io.py, orf.py, frameshift.py, report.py
   - output: report (file of type tbd)
  
   ## fasta.io Testing
   ```python
   from fasta_io import read_fasta

   records = read_fasta("example.fasta")
   for record in records:
     print(record["ID"], record["Sequence"])
   ```
