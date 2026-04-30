# ORF Frameshift Analysis & Comparison

## Team members

| Name               | GitHub username  | Email                        |
|--------------------|------------------|------------------------------|
| Emory Foerster     | emory-foerster   | efoerste@charlotte.edu       |
| Mekhi Lucas        | biomekhi         | mlucas37@charlotte.edu       |
| Akshitha Ganta     | akshitha-ganta   | aganta@charlotte.edu         |
| Runo Siakpebru     | runo313          | osiakpeb@charlotte.edu       |

## Description

This tool scans nucleotide FASTA sequences to identify Open Reading Frames (ORFs) and detect candidates that show likely frameshift events. For each sequence, it identifies all ORFs across reading frames, determines the longest ORF, and assesses whether a frameshift is present by detecting overlapping downstream ORFs in shifted frames.

The tool was designed and validated on coronavirus sequences (SARS-CoV-2, SARS-CoV-1) and currently supports gene annotation for coronaviruses. The core ORF detection and frameshift logic will work on any nucleotide FASTA input, but gene-level annotations in the HTML report are only available for coronavirus sequences.

The tool produces an interactive HTML report with five sections:

1. **Overview** — a summary dashboard showing total sequences analyzed, total ORFs found, frameshift count, and key biological findings with interpretation (e.g., identifying known programmed ribosomal frameshifts at the ORF1a/ORF1b boundary in coronaviruses).
2. **Per-sequence results** — for each input sequence: an ORF table with gene annotations, gene coverage visualization, frameshift position and shift type, dominance ratio, and sequence-specific context (e.g., SARS-CoV-2 reference genome annotations).
3. **Comparison** — a cross-sequence comparison showing how ORF counts and frameshift positions vary across all input sequences.
4. **Frameshift Plot** — an interactive Plotly visualization of frameshift candidates across sequences.
5. **Glossary** — definitions of key terms with literature references.

It also outputs a CSV summary table and a detailed JSON file with full ORF sequences and frameshift details per sequence.

## Biological relevance

Frameshifts are mutations that alter the reading frame of a gene, often producing truncated or non-functional proteins. Detecting potential frameshifts in nucleotide sequences is important for understanding viral evolution, identifying functional gene variants, and annotating novel sequences. This tool provides a lightweight, programmatic approach to flagging candidate frameshift sites without requiring alignment to a reference genome.

## Dependencies

- Python 3.10 or later
- pandas
- plotly
- colorema

The recommended way to set up the environment is with conda. To set up the environment and run the program enter the command below.

```bash
conda create -f environment.yml
```
Alternatively, you can install the core dependencies manually with:

```bash
pip install pandas plotly colorema 
```

## Usage

From the `src/` directory, run:

```bash
python3 main.py -f <fasta_file> [-v] [-m <min_orf_length>]
```

### Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-f`, `--fasta_file` | Yes | `../datasets/sequences.fasta` | Path to the input FASTA file. |
| `-r`, `--reference_fasta` | No | `../datasets/sequence_ref.fasta` | Path to the reference FASTA file for annotation. |
| `-m`, `--min_length` | No | `150` | Minimum ORF length in nucleotides. |
| `-oc`, `--output_csv` | No | — | File name for the output CSV file. |
| `-oj`, `--output_json` | No | — | File name for the output JSON file. |
| `-oh`, `--output_html` | No | — | File name for the output HTML report. |
| `-op`, `--output_plot` | No | — | File name for the frameshift track plot (HTML). |
| `-d`, `--output_dir` | No | `../examples` | Directory to save all output files. |
| `-v`, `--visualize` | No | — | Print results to terminal in addition to output files. |

### Example

```bash
python3 main.py -f ../datasets/sequences.fasta -v -m 150
```

## Output

The program produces the following output files in the `out_results/` directory:

| File | Description |
|------|-------------|
| `report.csv` | Table summarizing each sequence: frameshift status, position, shift type, dominance ratio, and longest ORF coordinates. |
| `report.json` | Detailed JSON output including full ORF sequences, frameshift details, and neighboring ORF information. |
| `full_report.html` | Interactive HTML report visualizing ORF structure and frameshift candidates per sequence. |
| `frameshift_plot.html` | Interactive HTML report visualizing frameshift plots per sequence. |
| `frame_job_out.txt` | Text file containing the script run information. |


## Frameshift detection logic

The tool detects potential frameshift events using the following approach:

1. **ORF detection** — all ORFs across all three reading frames are identified for each input sequence using a minimum length threshold (default: 150 nt).
2. **Longest ORF selection** — the longest ORF is selected as the primary candidate for frameshift analysis.
3. **Neighboring ORF search** — the tool searches within a 200 nucleotide window around the end position of the longest ORF for any ORFs that exist in a **different reading frame**. These are considered frameshift candidates.
4. **Shift type classification** — if neighboring ORFs are found, the frameshift type is calculated by subtracting the reading frames modulo 3, producing a `+1` or `+2` shift label.
5. **Dominance ratio** — the length of the longest ORF is expressed as a proportion of total nucleotides across all ORFs, giving a measure of how dominant it is in the sequence.

A frameshift is flagged as `True` if any neighboring ORF in a different frame is found within the search window. The position, frame, shift type, and full sequence of each neighboring ORF are recorded in the output.

## Limitations

- **Coronavirus-specific gene annotation** — gene-level annotations in the HTML report (e.g., Spike, ORF1a, N protein) are only available for coronavirus sequences. Other input sequences will still be analyzed for ORFs and frameshifts but will not receive gene labels.
- **False positive ORFs** — the tool identifies all sequences between start and stop codons above the minimum length threshold. In long sequences this produces many short, likely non-functional ORFs that may not represent real genes.
- **Frameshift detection sensitivity** — the neighboring ORF window approach may miss frameshifts that occur far from the end of the longest ORF, or produce false positives when unrelated ORFs happen to fall within the search window by chance.
- **Single longest ORF focus** — only the longest ORF per sequence is analyzed for frameshifts. Shorter ORFs that may contain biologically meaningful frameshifts are not evaluated.
- **Fixed window size** — the 200 nucleotide search window for neighboring ORFs is a fixed parameter. For sequences with very different genome sizes or gene densities this may need to be adjusted to avoid missing frameshifts or producing false positives.



A Bash test script is included to verify that the program runs correctly and produces the expected output. With the conda environment activated, run from the project root:

```bash
bash run_test.sh
```

The script runs `main.py` on a test dataset and checks that the expected output files are produced.

## Project structure

```
BINF_6112_Final/
├── README.md
├── LICENSE
├── ENVIRONMENT.md
├── dependencies.txt
├── run_test.sh
├── pseudocode/
│   ├── pseudocode.txt
│   └── flowcharts.txt
├── datasets/
│   ├── seq1.fasta
│   └── reference_sequence.fasta
├── src/
│   ├── __init__.py
│   ├── main.py
│   ├── fasta_io.py
│   ├── orf.py
│   ├── frameshift.py
│   ├── report.py
│   ├── visualize.py
│   └── html_report2.py
└── out_results/
    ├── report.csv
    ├── report.json
    ├── full_report.html
    ├── frameshift_plot.html
    └── frame_job_out.txt
```

## Module descriptions

| Module | Input | Output |
|--------|-------|--------|
| `fasta_io.py` | FASTA file | Dictionary of sequences |
| `orf.py` | Dictionary of sequences | List of longest ORFs per entry |
| `frameshift.py` | List of longest ORFs | Frameshift detection results |
| `report.py` | Frameshift results | CSV and JSON summary files |
| `visualize.py` / `html_report.py` | Report data | HTML visualization |
| `main.py` | All modules | Final output files |

## License

This project is licensed under the [GNU GPL V3 License](LICENSE). We chose this license to make our software accessible to all while ensuring that the source code and any modifications remain open under the same license.

## References

- Please see REFERENCES.MD
