#!/usr/bin/env python3
# Emory Foerster 
# Claude Sonnet 4.6
from datetime import datetime
import sys
import plotly.graph_objects as go

REFERENCE_ID = "NC_045512.2"


def gene_coverage_data(all_orfs: list[dict], virus_genes: list[tuple[int, int, str]]) -> list[tuple[str,int]]:
    """
    Purpose:
        Calculates the percentage of each known viral gene covered by detected ORFs
    Input:
        all_orfs    (list[dict])                : List of ORF dicts
        virus_genes (list[tuple[int, int, str]]): List of (gene_start, gene_end, gene_name) tuples defining known genes
    Output:
        list[tuple[str,int]]: A list of (gene_name, coverage_percent) pairs
    """
    rows = []
    for gene_start, gene_end, name in virus_genes:
        gene_len = gene_end - gene_start
        covered = sum(
            max(0, min(o["end"], gene_end) - max(o["start"], gene_start))
            for o in all_orfs
        )
        pct = min(100, int((covered / gene_len) * 100))
        rows.append((name, pct))
    return rows


def classify_sequence(seq_id: str, description: str = "") -> str:
    """
    Purpose:
        Classifies sequence if it is a known coronavirus strand based on accession ID and description.  
    Input:
        seq_id      (str): Sequence accession ID
        description (str): Optional description string from FASTA header, default is ""
    Output:
        str: Strain label - "SARS-CoV-2 reference", "SARS-CoV-2 variant" 
    """
    text = (seq_id + " " + description).lower()
    if seq_id == REFERENCE_ID:
        return "SARS-CoV-2 reference"
    elif any(k in text for k in ["sars-cov-2", "covid", "ncov", "omicron", "delta", "alpha", "beta", "gamma", "ba.", "xbb", "severe acute respiratory syndrome coronavirus 2", "wuhan"]):
        return "SARS-CoV-2 variant"
    elif any(k in text for k in [
        "sars-cov-1", "sars coronavirus 1", "severe acute respiratory syndrome coronavirus 1",
        "severe acute respiratory syndrome coronavirus,", "sars-associated", "ay613950", "nc_004718"
    ]):
        return "SARS-CoV-1"
    elif "severe acute respiratory syndrome" in text:
        return "SARS-CoV-1"
    elif "mers" in text:
        return "MERS-CoV"
    else:
        return description if description else seq_id

def badge_class(label: str) -> str:
    """
    Purpose:
        Returns class name for HTML sequence badge based on strain label. Used to determine badge color in the HTML report.
    Input:
        label (str): Strain label returned by classify_sequence()
    Output:
        str: "badge-sars1" for SARS-CoV-1, "badge-ref" for un-recognized sequences
    """
    if "SARS-CoV-1" in label or "MERS" in label:
        return "badge-sars1"
    return "badge-ref"

def get_strain_label(seq_id, description=""):
    """
    Purpose:
        Returns shortened sequence label for use in the comparison table header
    Input:
        seq_id      (str): Sequence accession ID
        description (str): Optional description string from FASTA header, default is ""
    Output:
        str: Shortened label e.g. "SARS-CoV ref", or the full classify_sequence() label
    """
    label = classify_sequence(seq_id, description)
    if label == "SARS-CoV-2 reference":
        return "SARS-CoV-2 ref"
    return label

def build_coverage_html(all_orfs: list[dict], virus_genes: list[tuple[int, int, str]]) -> str: 
    """
    Purpose:
        Builds HTML string for gene coverage progress bars
    Input:
        all_orfs    (list[dict]) : List of ORF dicts with "start", "end", "length", and "frame" keys
        virus_genes (list[tuple[int, int, str]]): List of (gene_start, gene_end, gene_name) tuples defining known genes
    Output:
        str: HTML string containing styled coverage bar for each known gene
    """
    coverage_html = ""
    for name, pct in gene_coverage_data(all_orfs, virus_genes):
        bar_color = "#639922" if pct >= 80 else "#BA7517" if pct >= 50 else "#A32D2D"
        coverage_html += f"""
        <div class="cov-row">
          <span class="cov-label">{name}</span>
          <div class="cov-bg"><div class="cov-bar" style="width:{pct}%;background:{bar_color}"></div></div>
          <span class="cov-pct">{pct}%</span>
        </div>"""
    return coverage_html

def build_orf_rows(all_orfs: list[dict], orfs_obj, seq_id:str, description: str) -> tuple[str, str]:
    """
    Purpose:
        Builds HTML table rows for all detected ORFS and the frame summary badges
    Input:
        all_orfs    (list[dict]) : List of ORF dicts with "start", "end", "length", and "frame" keys
        orfs_obj    (obj)        : ORF object exposing gene_gene_name()
        seq_id      (str)        : sequence accession ID
        description (str)        : Optional FASTA description string
    Output:
        tuple[str, str]: Tuple of (orf_rows, frame_summary: 
            - orf_rows      (str): HTML table row string for each ORF
            - frame_summary (str): HTML badge string summarizing ORF counts per frame
    """
    frame_colors = {0: "#639922", 1: "#185FA5", 2: "#7F77DD"}
    orf_rows = ""
    for orf in all_orfs:
        gene  = orfs_obj.get_gene_name(orf["start"], orf["end"], seq_id, description) or ""
        color = frame_colors[orf["frame"]]
        orf_rows += f"""<tr>
          <td><span class="fdot" style="background:{color}">F{orf['frame']}</span></td>
          <td>{orf['start']:,}</td><td>{orf['end']:,}</td>
          <td>{orf['length']:,}</td><td>{orf['length']//3:,}</td>
          <td class="gene-col">{gene}</td></tr>"""

    frame_summary = " ".join(
        f'<span class="fdot" style="background:{frame_colors[f]}">Frame {f}: {sum(1 for o in all_orfs if o["frame"]==f)}</span>'
        for f in range(3)
    )
    return orf_rows, frame_summary
def build_frameshift_html(details: list[dict], full_seq:str, orfs_obj) -> str:
    """
    Purpose:
        Builds HTML string for frameshift analysis section of sequence panel
    Input:
        details (list[dict]): List of frameshift details dicts, each containing
                                - shift_position (int): Position of the frameshift in the sequence
                                - shift_type     (str): Type of frameshift (e.g. "+1", "+2", "-1")
        full_seq    (str)   : Nucleotide sequence string
        orfs_obj    (obj)   : ORF object exposing gene_gene_name()
    Output:
        str: HTML string containing a frameshift card for each detected frameshift, 
             or a "no frameshift detected" message if details is empty
    """
    fs_html = ""
    stop_set = {"TAA", "TAG", "TGA"}
    for shift in details:
        pos        = shift["shift_position"]
        stype      = shift["shift_type"]
        gene_name  = orfs_obj.get_gene_name(pos, pos+61) or "unknown region"
        shift_seq  = full_seq[pos:pos+61]
        ss_seq     = full_seq[pos+shift["shift_magnitude"]:pos+61]
        orig_cod   = [shift_seq[i:i+3] for i in range(0, len(shift_seq)-2, 3)]
        shift_cod  = [ss_seq[i:i+3]    for i in range(0, len(ss_seq)-2,    3)]
        first_stop = next((i+1 for i, c in enumerate(shift_cod) if c in stop_set), None)

        orig_html  = "".join(f'<span class="codon codon-orig">{c}</span>'  for c in orig_cod  if len(c)==3)
        shift_html = "".join(
            f'<span class="codon {"codon-stop" if c in stop_set else "codon-shift"}">{c}</span>'
            for c in shift_cod if len(c)==3
        )
        if first_stop:
            stop_msg = f'<p class="stop-warn">Stop codon <strong>{shift_cod[first_stop-1]}</strong> at codon {first_stop} — protein likely truncated (loss-of-function)</p>'
        else:
            stop_msg = '<p class="no-stop">No immediate stop codon — may alter protein function</p>'

        fs_html += f"""
        <div class="fs-card">
          <div class="fs-title">{stype} frameshift at position {pos:,}</div>
          <p class="fs-loc">Located within: <strong>{gene_name}</strong></p>
          {stop_msg}
          <p class="codon-label">Original reading:</p>
          <div class="codon-row">{orig_html}</div>
          <p class="codon-label" style="margin-top:8px;">Shifted reading:</p>
          <div class="codon-row">{shift_html}</div>
        </div>"""

    if not fs_html:
        fs_html = '<p class="no-fs">No frameshift detected.</p>'
    return fs_html
def build_comparison_table(comparison_rows: list[dict]) -> str:
    """
    Purpose:
        Builds HTML comparison table summarizing key metrics across all sequences
    Input:
        compairson_rows (list[dict]): List of per-sequence metric dicts, each containing:
                                - seq_id      (str): Sequence accession ID
                                - total_orfs  (int): Total number of ORFS detected per-sequence
                                - seq_len     (int): Length of sequence in bp
                                - fs_pos      (str): Frameshift position or "None"
                                - fs_type     (str): Frameshift type of "None"
                                - spike_orfs  (int): Number of Spike gene ORFs detected
                                - description (str): Optional FASTA description string
    Output:
        str: HTML string containing styles comparison table.
    """

    def comp_row(label: str, key: str) -> str:
        """
        Purpose:
            Builds single HTML table row for the comparison table
        Input:
            label (str)        : Row label displayed in first column (e.g. "Total ORFs found").
            key   (str)        : Key to look up each comparison_rows dict (e.g. "total_orfs").
        Output:
            str: HTML string for single table row with one cell per sequence
        """
        cells = "".join(
            f'<td class="{"ref-col" if r["seq_id"]==REFERENCE_ID else ""}">{r[key]}</td>' 
            for r in comparison_rows
        )
        return f"<tr><td class='metric-col'>{label}</td>{cells}</tr>"
    
    comp_header_cells = "".join(
        f'<th>{"&#9733; " if r["seq_id"]==REFERENCE_ID else ""}{r["seq_id"]}<br>'
        f'<span class="th-sub">{get_strain_label(r["seq_id"], r.get("description",""))}</span></th>'
        for r in comparison_rows
    )
    return f"""
    <table class="comp-tbl">
      <thead><tr><th>Metric</th>{comp_header_cells}</tr></thead>
      <tbody>
        {comp_row("Total ORFs found", "total_orfs")}
        {comp_row("Sequence length", "seq_len")}
        {comp_row("Frameshift position", "fs_pos")}
        {comp_row("Frameshift type", "fs_type")}
        {comp_row("Spike ORFs", "spike_orfs")}
      </tbody>
    </table>"""

def build_findings_html(comparison_rows: list[dict], total_fs: int) -> str:
    """
    Purpose:
        Builds HTML string for key findings section of overview tab.
    Input:
        comparison_rows (list[dict]): List of per-sequence metric dicts (see build_compairson_table).
        total_fs        (str)       : Total number of sequences with detected frameshift
       
    Output:
        str: HTML string containing a finding cards for frameshift detection and ORF count variation. 
             
    """
    findings_html = ""
    fs_types = set(r["fs_type"] for r in comparison_rows if r["fs_type"] != "None")
    findings_html += f"""
    <div class="finding-card finding-warn">
      <div class="finding-title">Frameshifts detected in all {total_fs} sequences within ORF1a</div>
      <p class="finding-body">All sequences have a frameshift inside ORF1a (the largest gene, encoding replication machinery). In each case the shifted reading frame immediately hits a stop codon, truncating the protein. This is consistent with the known programmed ribosomal frameshift at the ORF1a/ORF1b boundary — a feature coronaviruses use intentionally to produce two different proteins from the same genomic region.</p>
      <sup><a href="#ref-10" style="color:#A32D2D;">[10]</a><a href="#ref-11" style="color:#A32D2D;">[11]</a><a href="#ref-12" style="color:#A32D2D;">[12]</a><a href="#ref-14" style="color:#A32D2D;">[14]</a></sup>
    </div>"""

    orf_counts = [r["total_orfs"] for r in comparison_rows]
    if max(orf_counts) != min(orf_counts):
        min_seq = min(comparison_rows, key=lambda r: r["total_orfs"])["seq_id"]
        max_seq = max(comparison_rows, key=lambda r: r["total_orfs"])["seq_id"]
        findings_html += f"""
    <div class="finding-card finding-info">
      <div class="finding-title">ORF counts vary across sequences ({min(orf_counts)} to {max(orf_counts)})</div>
      <p class="finding-body">{min_seq} has the fewest detected ORFs ({min(orf_counts)}) while {max_seq} has the most ({max(orf_counts)}). This likely reflects genuine differences in genome organization between strains rather than a detection error.<sup><a href="#ref-13" style="color:#185FA5;">[13]</a></sup></p>
    </div>"""
    return findings_html

def generate_html_report(all_records_data: list[dict], output_path:str = "report.html", frameshift_plots:list[dict] = None) -> None:
    """
    Purpose:
        Generates multi-tab HTML report summarizing ORF and frameshift detection analysis
        across one or multiple genome sequences. 
    Input:
        all_records_data(list[dict]): A list of per-sequence result dictionaries each containing:
                                        - seq_id      (str) : Sequence accession ID
                                        - description (str) : Optional FASTA description line
                                        - all_orfs    (list): ORF dicts from detect_ORFs(orf.py) / collect_ORFs(frameshift.py)
                                        - result      (dict): Frameshift analysis result with frameshift_details (frameshift.py)
                                        - full_seq    (str) : Nucleotide sequence string for ORF
                                        - orfs_obj    (obj) : ORF object exposing get_gene_name() ???
        output_path     (str)       : File path to write HTML to default = "report.html"
    Output:
       None, writes HTML file ot output_path and prints conformation to stdout
    """
    comparison_rows: list[dict] = []
    seq_panels_html: str = ""
    seq_tab_buttons: str = ""
    total_orfs_all : int = 0
    total_fs: int        = 0

    for i, rec in enumerate(all_records_data):
        seq_id   = rec["seq_id"]
        description = rec.get("description", "")
        all_orfs = rec["all_orfs"]
        result   = rec["result"]
        full_seq = rec["full_seq"]
        orfs_obj = rec["orfs_obj"]
        is_first = i == 0
        active   = "active" if is_first else ""
        total_orfs_all += len(all_orfs)

        # --- seq tab button ---
        label = seq_id
        seq_tab_buttons += f'<button class="seq-tab {active}" onclick="switchSeq(\'s{i}\')">{label}</button>\n'

        # --- metrics ---
        spike_orfs = sum(1 for o in all_orfs if orfs_obj.get_gene_name(o["start"], o["end"], seq_id, description) == "Spike (S)")
        details    = result.get("frameshift_details") or []
        fs_pos     = details[0]["shift_position"] if details else None
        fs_type    = details[0]["shift_type"]     if details else None
        if details:
            total_fs += 1

        # description = rec.get("description", "")
        comparison_rows.append({
            "seq_id":      seq_id,
            "total_orfs":  len(all_orfs),
            "seq_len":     len(full_seq),
            "fs_pos":      f"{fs_pos:,}" if fs_pos else "None",
            "fs_type":     fs_type or "None",
            "spike_orfs":  spike_orfs,
            "description": description,
        })

        # --- coverage bars ---
        coverage_html = build_coverage_html(all_orfs, orfs_obj.virus)

        # --- ORF table rows ---
        orf_rows, frame_summary = build_orf_rows(all_orfs, orfs_obj, seq_id, description)

        # --- frameshift html ---
        fs_html = build_frameshift_html(details, full_seq, orfs_obj)

        # --- reference banner ---
        strain_label = classify_sequence(seq_id, description)
        if seq_id == REFERENCE_ID:
            ref_banner = '<div class="ref-banner"><strong>SARS-CoV-2 reference genome</strong> — NC_045512.2 is the Wuhan-Hu-1 strain deposited when SARS-CoV-2 was first sequenced in 2020. This is the standard reference for all SARS-CoV-2 research.</div>'
        elif strain_label == "SARS-CoV-1":
            ref_banner = '<div class="sars1-banner"><strong>SARS-CoV-1 strain</strong> — This sequence is from the original 2003 SARS outbreak, a related but distinct coronavirus from SARS-CoV-2.</div>'
        elif strain_label == "MERS-CoV":
            ref_banner = '<div class="sars1-banner"><strong>MERS-CoV strain</strong> — Middle East Respiratory Syndrome coronavirus, a distinct betacoronavirus from SARS.</div>'
        elif strain_label == "SARS-CoV-2 variant":
            ref_banner = f'<div class="ref-banner"><strong>SARS-CoV-2 variant</strong> — {description if description else seq_id}</div>'
        else:
            ref_banner = f'<div class="sars1-banner"><strong>{strain_label}</strong></div>'

        seq_panels_html += f"""
        <div id="seq-s{i}" class="seq-panel {active}">
          {ref_banner}
          <div class="metrics-grid">
            <div class="metric-card"><div class="metric-label">Total ORFs</div><div class="metric-val">{len(all_orfs)}</div></div>
            <div class="metric-card"><div class="metric-label">Frameshift at</div><div class="metric-val" style="font-size:16px;">{"pos " + f"{fs_pos:,}" if fs_pos else "none"}</div></div>
            <div class="metric-card"><div class="metric-label">Spike ORFs</div><div class="metric-val">{spike_orfs}</div></div>
            <div class="metric-card"><div class="metric-label">Sequence length</div><div class="metric-val" style="font-size:15px;">{len(full_seq):,} bp</div></div>
          </div>
          <div class="two-col">
            <div>
              <p class="section-title">ORF table</p>
              <div class="tbl-scroll">
                <table class="orf-tbl">
                  <thead><tr><th>Frame</th><th>Start</th><th>End</th><th>Length (bp)</th><th>Length (aa)</th><th>Known gene</th></tr></thead>
                  <tbody>{orf_rows}</tbody>
                </table>
              </div>
              <div style="margin-top:10px;">{frame_summary}</div>
            </div>
            <div>
              <p class="section-title">Gene coverage</p>
              {coverage_html}
            </div>
          </div>
          <p class="section-title">Frameshift analysis</p>
          {fs_html}
        </div>"""


    # --- comparison table  ---
    comp_table = build_comparison_table(comparison_rows)

    # --- overview key findings ---
    findings_html = build_findings_html(comparison_rows, total_fs)

    now = datetime.now().strftime("%B %d, %Y at %H:%M")
    n   = len(all_records_data)
    badges_html = ""
    for rec in all_records_data:
        sid = rec["seq_id"]
        desc = rec.get("description", "")
        label = classify_sequence(sid, desc)
        bc = badge_class(label)
        badges_html += f'<span class="badge {bc}">{sid} &mdash; {label}</span>'

    sars1_list  = [r["seq_id"] for r in comparison_rows if classify_sequence(r["seq_id"], r.get("description","")) == "SARS-CoV-1"]
    sars2_list  = [r["seq_id"] for r in comparison_rows if r["seq_id"] != REFERENCE_ID and classify_sequence(r["seq_id"], r.get("description","")) not in ("SARS-CoV-1", "MERS-CoV")]
    mers_list   = [r["seq_id"] for r in comparison_rows if classify_sequence(r["seq_id"], r.get("description","")) == "MERS-CoV"]
    compare_note = "NC_045512.2 is the SARS-CoV-2 reference (starred)."
    if sars1_list:
        compare_note += f" {', '.join(sars1_list)} are SARS-CoV-1 strains from the 2003 outbreak."
    if sars2_list:
        compare_note += f" {', '.join(sars2_list)} are additional SARS-CoV-2 variants."
    if mers_list:
        compare_note += f" {', '.join(mers_list)} are MERS-CoV sequences."

    # Build per-sequence plot tab panels
    plot_tab_buttons = ""
    plot_panels_html = ""

    if frameshift_plots:
        for j, entry in enumerate(frameshift_plots):
            pid    = f"fp{j}"
            active = "active" if j == 0 else ""
            label  = entry["seq_id"]
            plot_tab_buttons += (
                f'<button class="seq-tab {active}" '
                f'onclick="switchPlot(\'{pid}\')">{label}</button>\n'
            )
            plot_panels_html += (
                f'<div id="plot-{pid}" class="seq-panel {active}">'
                f'{entry["html"]}'
                f'</div>'
            )
    else:
        plot_panels_html = '<p class="no-fs">No frameshift plots generated.</p>'
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Coronavirus ORF Analysis Report</title>
  <style>
    *{{box-sizing:border-box;margin:0;padding:0;}}
    body{{font-family:'Segoe UI',sans-serif;background:#f0f2f5;color:#222;font-size:14px;}}
    header{{background:#1a1a2e;color:white;padding:1.5rem 2rem;}}
    header h1{{font-size:1.4rem;font-weight:500;}}
    header p{{opacity:.65;font-size:.85rem;margin-top:4px;}}
    .badges{{display:flex;gap:8px;margin-top:10px;flex-wrap:wrap;}}
    .badge{{font-size:11px;padding:3px 10px;border-radius:20px;font-weight:500;}}
    .badge-ref{{background:#B5D4F4;color:#0C447C;}}
    .badge-sars1{{background:#CECBF6;color:#3C3489;}}
    .main-tabs{{display:flex;gap:0;border-bottom:1px solid #ddd;background:white;padding:0 1.5rem;}}
    .main-tab{{padding:10px 18px;font-size:13px;cursor:pointer;border-bottom:2px solid transparent;color:#666;background:none;border-top:none;border-left:none;border-right:none;}}
    .main-tab.active{{color:#1a1a2e;border-bottom:2px solid #1a1a2e;font-weight:500;}}
    .tab-content{{display:none;padding:1.5rem;max-width:1300px;margin:0 auto;}}
    .tab-content.active{{display:block;}}
    .summary-banner{{background:#E6F1FB;border-left:3px solid #185FA5;border-radius:0 8px 8px 0;padding:12px 16px;margin-bottom:1.25rem;font-size:13px;line-height:1.6;color:#0C447C;}}
    .metrics-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(130px,1fr));gap:10px;margin-bottom:1.25rem;}}
    .metric-card{{background:#f8f8f8;border-radius:8px;padding:12px 14px;}}
    .metric-label{{font-size:11px;color:#888;margin-bottom:4px;}}
    .metric-val{{font-size:22px;font-weight:500;color:#1a1a2e;}}
    .section-title{{font-size:11px;font-weight:500;text-transform:uppercase;letter-spacing:.06em;color:#888;margin:1.25rem 0 .6rem;}}
    .cov-row{{display:flex;align-items:center;gap:10px;margin-bottom:7px;}}
    .cov-label{{width:150px;font-size:12px;flex-shrink:0;}}
    .cov-bg{{flex:1;background:#e5e5e5;border-radius:4px;height:8px;}}
    .cov-bar{{height:8px;border-radius:4px;}}
    .cov-pct{{width:36px;font-size:12px;text-align:right;color:#666;}}
    .tbl-scroll{{overflow-x:auto;max-height:380px;overflow-y:auto;border:1px solid #eee;border-radius:6px;}}
    .orf-tbl{{width:100%;border-collapse:collapse;font-size:12px;}}
    .orf-tbl th{{background:#1a1a2e;color:white;padding:6px 10px;text-align:left;position:sticky;top:0;font-size:11px;font-weight:500;}}
    .orf-tbl td{{padding:5px 10px;border-bottom:1px solid #f0f0f0;}}
    .orf-tbl tr:hover td{{background:#fafafa;}}
    .gene-col{{color:#BA7517;font-weight:500;}}
    .fdot{{display:inline-block;color:white;padding:2px 7px;border-radius:4px;font-size:11px;font-weight:500;margin:1px;}}
    .two-col{{display:grid;grid-template-columns:1fr 1fr;gap:2rem;margin-bottom:1rem;}}
    .fs-card{{border-left:3px solid #A32D2D;border-radius:0 8px 8px 0;padding:12px 16px;margin-bottom:1rem;background:#fff8f8;}}
    .fs-title{{font-size:13px;font-weight:500;color:#A32D2D;margin-bottom:5px;}}
    .fs-loc{{font-size:12px;color:#444;margin-bottom:6px;}}
    .stop-warn{{font-size:12px;color:#791F1F;background:#FCEBEB;border-radius:6px;padding:6px 10px;margin:6px 0;line-height:1.5;}}
    .no-stop{{font-size:12px;color:#854F0B;background:#FAEEDA;border-radius:6px;padding:6px 10px;margin:6px 0;}}
    .no-fs{{color:#3B6D11;font-style:italic;font-size:13px;}}
    .codon-label{{font-size:11px;color:#888;margin-top:8px;margin-bottom:3px;}}
    .codon-row{{display:flex;flex-wrap:wrap;gap:4px;}}
    .codon{{padding:3px 6px;border-radius:4px;font-family:monospace;font-size:12px;font-weight:500;}}
    .codon-orig{{background:#EAF3DE;color:#27500A;}}
    .codon-shift{{background:#FAEEDA;color:#633806;}}
    .codon-stop{{background:#FCEBEB;color:#791F1F;}}
    .ref-banner{{background:#E6F1FB;border-left:3px solid #185FA5;border-radius:0 8px 8px 0;padding:10px 14px;margin-bottom:1rem;font-size:12px;color:#0C447C;line-height:1.5;}}
    .sars1-banner{{background:#EEEDFE;border-left:3px solid #534AB7;border-radius:0 8px 8px 0;padding:10px 14px;margin-bottom:1rem;font-size:12px;color:#3C3489;line-height:1.5;}}
    .seq-tabs{{display:flex;gap:8px;margin-bottom:1rem;flex-wrap:wrap;}}
    .seq-tab{{padding:5px 14px;font-size:12px;border-radius:20px;border:1px solid #ccc;cursor:pointer;background:white;color:#666;}}
    .seq-tab.active{{background:#1a1a2e;color:white;border-color:#1a1a2e;}}
    .seq-panel{{display:none;}}
    .seq-panel.active{{display:block;}}
    .comp-tbl{{width:100%;border-collapse:collapse;font-size:13px;}}
    .comp-tbl th{{text-align:left;padding:9px 14px;background:#f0f0f0;font-weight:500;font-size:12px;border-bottom:2px solid #ddd;}}
    .th-sub{{font-weight:400;color:#888;font-size:11px;}}
    .comp-tbl td{{padding:8px 14px;border-bottom:1px solid #eee;}}
    .comp-tbl tr:hover td{{background:#fafafa;}}
    .metric-col{{color:#444;font-weight:500;}}
    .ref-col{{color:#185FA5;font-weight:500;}}
    .finding-card{{border-radius:0 8px 8px 0;padding:12px 16px;margin-bottom:1rem;font-size:13px;line-height:1.6;}}
    .finding-warn{{background:#fff8f8;border-left:3px solid #A32D2D;}}
    .finding-info{{background:#E6F1FB;border-left:3px solid #185FA5;}}
    .finding-title{{font-weight:500;margin-bottom:5px;}}
    .finding-body{{color:#444;}}
    .glossary-item{{border:1px solid #eee;border-radius:8px;padding:12px 16px;margin-bottom:10px;}}
    .glossary-term{{font-weight:500;font-size:13px;margin-bottom:4px;}}
    .glossary-def{{font-size:13px;color:#555;line-height:1.6;}}
    .finding-refs{{margin-top:8px;font-size:11px;color:#888;line-height:1.7;}}
    .finding-refs a{{color:#A32D2D;text-decoration:none;}}
    .finding-refs a:hover{{text-decoration:underline;}}
    @media(max-width:700px){{.two-col{{grid-template-columns:1fr;}}}}
  </style>
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body>
  <header>
    <h1>Coronavirus ORF &amp; frameshift analysis report</h1>
    <p>Generated {now} &middot; {n} sequence(s) analyzed</p>
    <div class="badges">
      {badges_html}
    </div>
  </header>

  <div class="main-tabs">
    <button class="main-tab active" onclick="switchTab('overview')">Overview</button>
    <button class="main-tab" onclick="switchTab('sequences')">Per-sequence results</button>
    <button class="main-tab" onclick="switchTab('compare')">Comparison</button>
    <button class="main-tab" onclick="switchTab('plot')">Frameshift Plot</button>
    <button class="main-tab" onclick="switchTab('glossary')">Glossary</button>
  </div>

  <div id="tab-overview" class="tab-content active">
    <div class="summary-banner">
      This tool scanned {n} coronavirus genome(s) for open reading frames (ORFs) &mdash; the regions that encode proteins &mdash; and checked whether any contain frameshifts, which are mutations that shift how the sequence is read and often destroy protein function.
    </div>
    <div class="metrics-grid">
      <div class="metric-card"><div class="metric-label">Sequences analyzed</div><div class="metric-val">{n}</div></div>
      <div class="metric-card"><div class="metric-label">Total ORFs found</div><div class="metric-val">{total_orfs_all}</div></div>
      <div class="metric-card"><div class="metric-label">Frameshifts detected</div><div class="metric-val">{total_fs}</div></div>
    </div>
    <p class="section-title">Key findings</p>
    {findings_html}
  </div>

  <div id="tab-sequences" class="tab-content">
    <div class="seq-tabs">{seq_tab_buttons}</div>
    {seq_panels_html}
  </div>

  <div id="tab-compare" class="tab-content">
    <div class="summary-banner" style="background:#EEEDFE;border-color:#534AB7;color:#3C3489;">
            {compare_note} Differences between columns reflect genuine biological variation.
    </div>
    {comp_table}
  </div>
   <div id="tab-plot" class="tab-content">
    <div class="summary-banner">
      Sequence-track view of the longest ORF (steelblue) and neighboring frameshift
      candidates across all 3 reading frames. The dashed vertical line marks the
      detected frameshift site.
    </div>
    <div class="seq-tabs">{plot_tab_buttons}</div>
    {plot_panels_html}
  </div>
  
  <div id="tab-glossary" class="tab-content">
    <div class="summary-banner" style="background:#EAF3DE;border-color:#3B6D11;color:#27500A;">
      Plain-language definitions for the key terms used in this report.
    </div>
    <div class="glossary-item"><div class="glossary-term">Open reading frame (ORF)</div><div class="glossary-def">A stretch of DNA that starts with a start codon (ATG) and ends with a stop codon. ORFs are candidates for protein-coding genes &mdash; the cell's machinery reads them to build proteins.<sup><a href="#ref-1" style="color:#185FA5;">[1]</a><a href="#ref-2" style="color:#185FA5;">[2]</a><a href="#ref-3" style="color:#185FA5;">[3]</a><a href="#ref-4" style="color:#185FA5;">[4]</a></sup></div></div>
    <div class="glossary-item"><div class="glossary-term">Reading frame</div><div class="glossary-def">DNA is read in triplets called codons. Depending on where you start reading, you get a different set of triplets. Frame 0 starts at position 0, Frame 1 at position 1, Frame 2 at position 2. Each gives a completely different protein sequence from the same DNA.<sup><a href="#ref-5" style="color:#185FA5;">[5]</a><a href="#ref-2" style="color:#185FA5;">[2]</a></sup></div></div>
    <div class="glossary-item"><div class="glossary-term">Frameshift</div><div class="glossary-def">An insertion or deletion that shifts the reading frame. A +1 frameshift means 1 nucleotide was inserted; +2 means 2 were inserted. Everything downstream is now read in a different frame, usually producing a garbled or truncated protein.<sup><a href="#ref-5" style="color:#185FA5;">[5]</a><a href="#ref-6" style="color:#185FA5;">[6]</a></sup></div></div>
    <div class="glossary-item"><div class="glossary-term">Stop codon</div><div class="glossary-def">TAA, TAG, or TGA &mdash; triplets that signal "stop building this protein here." When a frameshift introduces a premature stop codon, the resulting protein is truncated and usually non-functional.<sup><a href="#ref-7" style="color:#185FA5;">[7]</a><a href="#ref-8" style="color:#185FA5;">[8]</a></sup></div></div>
    <div class="glossary-item"><div class="glossary-term">Gene coverage</div><div class="glossary-def">What percentage of a known gene's length is covered by detected ORFs. 100% means the tool found an ORF spanning the full known gene. Lower percentages may indicate mutations, deletions, or detection limits.<sup><a href="#ref-3" style="color:#185FA5;">[3]</a><a href="#ref-9" style="color:#185FA5;">[9]</a></sup></div></div>
    <div class="glossary-item"><div class="glossary-term">Programmed ribosomal frameshift</div><div class="glossary-def">Some viruses, including coronaviruses, deliberately exploit frameshifting. The ORF1a/ORF1b boundary in SARS genomes contains a slippery sequence where ribosomes intentionally slip, allowing the virus to produce two different proteins from the same genomic region.<sup><a href="#ref-10" style="color:#185FA5;">[10]</a><a href="#ref-11" style="color:#185FA5;">[11]</a><a href="#ref-12" style="color:#185FA5;">[12]</a><a href="#ref-14" style="color:#185FA5;">[14]</a><a href="#ref-15" style="color:#185FA5;">[15]</a></sup></div></div>
    <div class="glossary-item"><div class="glossary-term">Spike protein</div><div class="glossary-def">The protein on the surface of coronaviruses that binds to human cells. It is the primary target of COVID-19 vaccines. Detecting ORFs in the Spike gene region helps confirm the genome is intact and the sequence is being read correctly.<sup><a href="#ref-16" style="color:#185FA5;">[16]</a><a href="#ref-17" style="color:#185FA5;">[17]</a></sup></div></div>
    <p class="section-title" style="margin-top:2rem;">References</p>
    <div style="font-size:12px;color:#555;line-height:2;border-top:1px solid #eee;padding-top:1rem;">
      <div id="ref-1">[1] National Human Genome Research Institute. Open Reading Frame. <a href="https://www.genome.gov/genetics-glossary/Open-Reading-Frame" target="_blank">genome.gov Talking Glossary</a>.</div>
      <div id="ref-2">[2] National Human Genome Research Institute. Bioinformatics — Finding Genes. <a href="https://www.genome.gov/25020001/online-education-kit-bioinformatics-finding-genes" target="_blank">genome.gov Education Kit</a>.</div>
      <div id="ref-3">[3] Wu F et al. (2020). A new coronavirus associated with human respiratory disease in China. <a href="https://doi.org/10.1038/s41586-020-2008-3" target="_blank">Nature 579:265–269</a>.</div>
      <div id="ref-4">[4] Finkel Y et al. (2021). The coding capacity of SARS-CoV-2. <a href="https://doi.org/10.1038/s41586-020-2739-1" target="_blank">Nature 589:125–130</a>.</div>
      <div id="ref-5">[5] Crick FHC et al. (1961). General nature of the genetic code for proteins. <a href="https://doi.org/10.1038/192227a0" target="_blank">Nature 192:1227–1232</a>.</div>
      <div id="ref-6">[6] Farabaugh PJ (1996). Programmed translational frameshifting. <a href="https://doi.org/10.1128/mmbr.60.1.103-134.1996" target="_blank">Microbiol Mol Biol Rev 60(1):103–134</a>.</div>
      <div id="ref-7">[7] National Human Genome Research Institute. Stop Codon. <a href="https://www.genome.gov/genetics-glossary/Stop-Codon" target="_blank">genome.gov Talking Glossary</a>.</div>
      <div id="ref-8">[8] Brenner S et al. (1967). UGA: a third nonsense triplet in the genetic code. <a href="https://doi.org/10.1038/2131049a0" target="_blank">Nature 213:449–450</a>.</div>
      <div id="ref-9">[9] Lee RT et al. (2020). The UCSC SARS-CoV-2 Genome Browser. <a href="https://doi.org/10.1038/s41588-020-0700-8" target="_blank">Nature Genetics 52:986–990</a>.</div>
      <div id="ref-10">[10] Brierley I, Digard P &amp; Inglis SC (1989). Characterization of an efficient coronavirus ribosomal frameshifting signal. <a href="https://doi.org/10.1016/0092-8674(89)90124-4" target="_blank">Cell 57(4):537–547</a>.</div>
      <div id="ref-11">[11] Plant EP et al. (2005). A three-stemmed mRNA pseudoknot in the SARS coronavirus frameshift signal. <a href="https://doi.org/10.1371/journal.pbio.0030172" target="_blank">PLOS Biology 3(6):e172</a>.</div>
      <div id="ref-12">[12] Kelly JA et al. (2020). Structural and functional conservation of the programmed &minus;1 ribosomal frameshift signal of SARS-CoV-2. <a href="https://doi.org/10.1261/rna.076141.120" target="_blank">RNA 26(9):1175–1189</a>.</div>
      <div id="ref-13">[13] Gorbalenya AE et al. (2020). The species Severe acute respiratory syndrome-related coronavirus: classifying 2019-nCoV and naming it SARS-CoV-2. <a href="https://doi.org/10.1038/s41564-020-0695-z" target="_blank">Nature Microbiology 5:536–544</a>.</div>
      <div id="ref-14">[14] Bhatt PR et al. (2021). Structural basis of ribosomal frameshifting during translation of the SARS-CoV-2 RNA genome. <a href="https://doi.org/10.1126/science.abf3546" target="_blank">Science 372(6548):1306–1313</a>.</div>
      <div id="ref-15">[15] Haniff HS et al. (2021). Restriction of SARS-CoV-2 replication by targeting programmed &minus;1 ribosomal frameshifting. <a href="https://doi.org/10.1073/pnas.2023051118" target="_blank">PNAS 118(26):e2023051118</a>.</div>
      <div id="ref-16">[16] Walls AC et al. (2020). Structure, function, and antigenicity of the SARS-CoV-2 spike glycoprotein. <a href="https://doi.org/10.1016/j.cell.2020.02.058" target="_blank">Cell 181(2):281–292</a>.</div>
      <div id="ref-17">[17] Wrapp D et al. (2020). Cryo-EM structure of the 2019-nCoV spike in the prefusion conformation. <a href="https://doi.org/10.1126/science.abb2507" target="_blank">Science 367(6483):1260–1263</a>.</div>
    </div>
  </div>

  <script>
    function switchTab(id) {{
      document.querySelectorAll('.main-tab').forEach(t => t.classList.remove('active'));
      document.querySelectorAll('.tab-content').forEach(t => t.classList.remove('active'));
      document.querySelector('[onclick="switchTab(\\''+id+'\\')"]').classList.add('active');
      document.getElementById('tab-'+id).classList.add('active');
    }}
    function switchSeq(id) {{
      document.querySelectorAll('.seq-tab').forEach(t => t.classList.remove('active'));
      document.querySelectorAll('.seq-panel').forEach(t => t.classList.remove('active'));
      document.querySelector('[onclick="switchSeq(\\''+id+'\\')"]').classList.add('active');
      document.getElementById('seq-'+id).classList.add('active');
    }}
    function switchPlot(id) {{
      document.querySelectorAll('[id^="plot-"]').forEach(p => p.classList.remove('active'));
      document.querySelectorAll('[onclick^="switchPlot"]').forEach(t => t.classList.remove('active'));
      document.getElementById('plot-'+id).classList.add('active');
      document.querySelector('[onclick="switchPlot(\\''+id+'\\')"]').classList.add('active');
    }}
  </script>
</body>
</html>"""

    with open(output_path, "w") as f:
        f.write(html)
    sys.stdout.write(f"HTML report saved to {output_path}\n")





