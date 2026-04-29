import re
import sys
from typing import Optional
from colorama import Fore, Back, Style, init
from orf import parse_codons
 
init(autoreset=True)
 
# Known SARS-CoV-2 gene positions (start, end, name)
# Coordinates are 0-based and based on NC_045512.2 (Wuhan-Hu-1 reference)
sars_cov2_genes = [
    (266,   21555,  "ORF1a"),
    (266,   13468,  "nsp1-10"),
    (13468, 21555,  "nsp12-16 (RdRp)"),
    (21563, 25384,  "Spike (S)"),
    (25393, 26220,  "ORF3a"),
    (26245, 26472,  "E protein"),
    (26523, 27191,  "M protein"),
    (27202, 27387,  "ORF6"),
    (27394, 27759,  "ORF7a"),
    (27756, 27887,  "ORF7b"),
    (27894, 28259,  "ORF8"),
    (28274, 29533,  "N protein"),
    (29558, 29674,  "ORF10"),
]
class ORFS:
    """
    Purpose: 
        Visualization and annotation class for ORF detection results.
        Provides methods to annotate ORFs with known gene names
        Prints colorized terminal summaries, and display frameshift analysis.
    Attributes: 
        virus(list[tuple[int, int, str]]): List of (start, end, name) tuples represent known gene pos for reference virus
    """

    def __init__(self, virus_genes: list[tuple[int, int, str]]) -> None:
        """
        Purpose: 
            Initialize the ORFs visualizer with set of known gene coordinates
        Input: 
            virus_genes (list[tuple[int, int, str]]): List of (start, end, name) tuples with known gene positions
        Output: 
            None
        """
        self.virus = virus_genes
    def get_gene_name(self, orf_start: int, orf_end: int, seq_id: str = "", description: str="") -> Optional[str]:
        """
        Purpose: 
            Return name of known gene that best overlaps a given ORF
            ORF assigned to a gene if more than 50% of ORF's length overlaps with gene coordinate
        Inputs: 
            orf_start   (int): Start position of the ORF (0-based)
            orf_end     (int): End position of the ORF (0-based, exclusive)
            seq_id      (str): Sequence identifier
            description (str): FASTA header description 
        Output: 
            str or None: the name of the overlapping gene, or None if no gene overlaps by more than 50%
        """
        text = (seq_id + " " + description).lower()

        #known non-coronavirus viruses - skip annotation
        known_other = any(k in text for k in ["west nile", "influenze", "hiv", "dengue", "zika", "ebola", "rabies", "hepatitis"])
        if known_other: 
            return None

        start = orf_start
        end = orf_end
        length = end-start
        for gene_start, gene_end, name in sars_cov2_genes:
            # overlap starts at the start that is greater
            # overlap ends at the end that is smaller
            overlap_len = max(0, min(end, gene_end) - max(start,gene_start))
            # 0.5 = at least 50% of this ORF must overlap with known gene for us to label it as gene
            if overlap_len > 0 and overlap_len / length > 0.5:
                return name
        return None

    def visualize_orf(self, all_orfs: list[dict], seq_id: str, seq:str, description: str = "", min_length: int = 150) -> None:
        """
        Purpose: 
            Print a formatted, colorized ORF table and reading frame summary to terminal
            Each ORF is color-coded by frame and annotated with known gene if it exists
        Input: 
            all_orfs (list[dict]): List of ORF dicts, each containing keys: 
                                    'frame' (int), 'start'(int), 'end'(int), 'length'(int)
            seq_id      (str):     Sequence identifier (e.g "NC_045512.2")
            seq         (str):     Full nucleotide sequence string
            description (str):     FASTA header description string, default is "" 
        Output:
            None
        """
        frame_colors = {0: Back.GREEN, 1: Back.BLUE, 2: Back.MAGENTA}
        frame_text_colors = {0: Fore.BLACK, 1: Fore.WHITE, 2: Fore.BLACK}
        sys.stdout.write(f"\n{'='*70}\n")
        sys.stdout.write(f"  SARS-CoV-2 ORF Analysis — All 3 Reading Frames\n")
        sys.stdout.write(f"  Sequence length: {len(seq):,} bp\n")
        sys.stdout.write(f"  Sequence ID: {seq_id}\n")
        sys.stdout.write(f"  Description: {description}\n")
        sys.stdout.write(f"{'='*70}\n\n")

        sys.stdout.write(f"Found {len(all_orfs)} ORFs (≥{min_length}bp) across all 3 frames:\n\n")

        sys.stdout.write(f"{'Frame':<8}{'Start':<10}{'End':<10}{'Length (bp)':<14}{'Length (aa)':<14}{'Known Gene'}\n")
        sys.stdout.write("-" * 75 + "\n")
     
        for orf in all_orfs:
            gene = self.get_gene_name(orf["start"], orf["end"], seq_id, description)
            gene_str = gene if gene else ""
            color = frame_colors[orf["frame"]]
            text_color = frame_text_colors[orf['frame']]
            sys.stdout.write(
                f"{color}{text_color} F{orf['frame']} {Style.RESET_ALL}"
                f"  {orf['start']:<10}{orf['end']:<10}{orf['length']:<14}"
                f"{orf['length']//3:<14}"
                f"{Fore.YELLOW + gene_str + Style.RESET_ALL if gene_str else ''}\n"
            )
     
        # Summary by frame
        sys.stdout.write(f"\n{'='*70}\n")
        sys.stdout.write("  Summary by Reading Frame:\n")
        sys.stdout.write(f"{'='*70}\n")
        for frame in range(0, 3):
            frame_orfs = [o for o in all_orfs if o["frame"] == frame]
            color = frame_colors[frame]
            text_color = frame_text_colors[frame]

            sys.stdout.write(f"  {color}{text_color} Frame {frame} {Style.RESET_ALL}: {len(frame_orfs)} ORFs found\n")

    def display_frameshift(self, result: dict, seq:str, seq_id: str = "", description: str = "") -> None:
        """
        Purpose:
            Print biological interpretation of each detected frameshift to terminal 
            Includes: affected gene, if premature stop codon introduced, and colorized original vs shifted codon seq
        Input: 
            result(dict): Analysis result dict from FrameshiftDetector.analyze(), 
                          Must contain key 'frameshift_details' a list of shift dicts with keys: 
                            - 'shift_position' (int): Genomic position of shift
                            - 'shift_type'     (str): e.g "+1" or "+2"
                            -  shift_magnitude (int): Nucleotide shifted
            seq        (str): Full nucleotide seq string
            seq_id     (str): Sequence identifier, default is ""
            descriptio (str): FASTA header description string, default is "" 
        Output: 
            None
        """
        details = result["frameshift_details"]
        for shift in details:
            pos = shift['shift_position']
            shift_type = shift['shift_type']
            shift_seq = seq[pos:pos+61]
            ss_seq = seq[pos+shift["shift_magnitude"]: pos+61]
            codons = parse_codons(shift_seq, 0)
            shifted_codons = parse_codons(ss_seq, 0)

            # --- biological context ---
            gene_name = self.get_gene_name(pos, pos + 61, seq_id, description) or "unknown region"
            stop_codons = {"TAA", "TAG", "TGA"}
            shifted_codon_seqs = [item[0] for item in shifted_codons]
            first_stop = next(
                (i+1 for i, c in enumerate(shifted_codon_seqs) if c in stop_codons),
                None
            )

            sys.stdout.write(f"\n{'='*70}\n")
            sys.stdout.write(f"  *  Frameshift detected at position {pos} -- type {shift_type}\n")
            sys.stdout.write(f"  → Located within: {gene_name}\n")
            if first_stop:
                stop_codon = shifted_codon_seqs[first_stop - 1]
                sys.stdout.write(f"  → Shifted frame introduces STOP codon ({stop_codon}) at codon {first_stop}\n")
                sys.stdout.write(f"  → Likely truncates the protein — possible loss-of-function mutation\n")
            else:
                sys.stdout.write(f"  → No immediate stop codon in shifted frame (may alter protein function)\n")
            sys.stdout.write(f"{'='*70}\n")

            # ... rest of your existing codon printing code ...
            sys.stdout.write(f"\n  Original protein start (position {pos}):\n")
            for codon in [item[0] for item in codons]:
                sys.stdout.write(Back.GREEN + codon + Style.RESET_ALL + " ")
            sys.stdout.write("\n\n")

            sys.stdout.write(f"\n  Shifted protein start (position {pos}):\n")
            for i, codon in enumerate(shifted_codon_seqs):
                color = Back.RED if codon in stop_codons else Back.YELLOW
                sys.stdout.write(color + Fore.WHITE + codon + Style.RESET_ALL + " ")
            sys.stdout.write("\n\n")
 

    def display_gene_coverage(self, all_orfs: list[dict]) -> None:
        """
        Purpose: 
            Print colorized gene coverage bar chart to the terminal showing what percentage of each known gene it covered by detected ORFs
        Input: 
            all_orfs (list[dict]): List of ORF dicts, each contianing keys: 'start' (int) and 'end' (int)
        Output: 
            None - prints to stdout
                   Green bar   = > 80% covered
                   Yellow bar  = > 50% covered
                   Red bar     = < 50% covered
        """
        sys.stdout.write(f"\n{'='*70}\n")
        sys.stdout.write("  Gene Coverage Summary:\n")
        sys.stdout.write(f"{'='*70}\n")
        
        for gene_start, gene_end, name in self.virus:
            gene_len = gene_end - gene_start
            covered = 0
            for orf in all_orfs:
                overlap = max(0, min(orf["end"], gene_end) - max(orf["start"], gene_start))
                covered += overlap
            # cap at 100%
            pct = min(100, int((covered / gene_len) * 100))
            bar_filled = int(pct / 10)
            bar = "█" * bar_filled + "░" * (10 - bar_filled)
            color = Fore.GREEN if pct >= 80 else Fore.YELLOW if pct >= 50 else Fore.RED
            sys.stdout.write(f"  {name:<25} {color}{bar}{Style.RESET_ALL}  {pct}% covered \n")

def display_cross_sequence_comparison(all_results: list[dict]) -> None:
    """
    Purpose: 
        Print a formatted cross-sequence comparison table to terminal. 
        Compares key metrics across all analyzed sequences side by side.
    Input:
        all_results (list[dict]): List of result dicts, one per sequence,
                                  each containing:
                                      - 'sequence_id'    (str):      Accession ID.
                                      - 'total_orfs'     (int):      Total ORFs found.
                                      - 'frameshift_pos' (int|str):  Frameshift position
                                                                      or "None".
                                      - 'frameshift_type'(str):      e.g. "+1", "+2",
                                                                      or "None".
                                      - 'spike_orfs'     (int):      ORFs overlapping
                                                                      the Spike gene.
    Output: 
        None
    """
    sys.stdout.write(f"\n{'='*70}\n")
    sys.stdout.write("  Cross-Sequence Comparison (SARS-CoV-1 vs SARS-CoV-2):\n")
    sys.stdout.write(f"  Reference: NC_045512.2 (SARS-CoV-2 Wuhan-Hu-1)\n")
    sars1_ids = [r["sequence_id"] for r in all_results if r["sequence_id"] != "NC_045512.2"]

    sys.stdout.write(f"  Outgroups:  {', '.join(sars1_ids)}\n")
    sys.stdout.write(f"{'='*70}\n")

    ids    = [r["sequence_id"]    for r in all_results]
    col_w  = max(len(i) for i in ids) + 2

    header = f"  {'Metric':<30}" + "".join(f"{i:>{col_w}}" for i in ids)
    sys.stdout.write(header + "\n")
    sys.stdout.write("-" * len(header) + "\n")

    metrics = [
        ("Total ORFs found",   "total_orfs"),
        ("Frameshift position","frameshift_pos"),
        ("Frameshift type",    "frameshift_type"),
        ("Spike ORFs found",   "spike_orfs"),
    ]
    for label, key in metrics:
        row = f"  {label:<30}" + "".join(f"{str(r.get(key,'N/A')):>{col_w}}" for r in all_results)
        sys.stdout.write(row + "\n")
    sys.stdout.write(f"{'='*70}\n\n")

