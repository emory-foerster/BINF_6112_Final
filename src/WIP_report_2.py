#!/usr/bin/env python3
# Mekhi Lucas
from typing import List, Dict, Any, Optional
import pandas as pd
import os
import sys
import argparse
import plotly.graph_objects as go
from fasta_io import read_fasta
from orf import detect_all_frames
from frameshift import FrameshiftDetector


class OrfReport:
    """
    Generate output files (JSON, CSV, HTML, frameshift plot) from ORF + frameshift results.

    Purpose:
        Take analysis results and write them to one or more report formats.

    Input:
        results   (List[Dict[str, Any]]): list of per-sequence result dictionaries.
        output_dir (str):                 directory where report files will be saved.
        all_orfs  (List[Dict] | None):    full ORF list from orf.detect_all_frames —
                                          required for neighbor dominance in the
                                          frameshift plot; safe to omit for other formats.

    Output:
        Report files produced in output_dir.
    """

    def __init__(
        self,
        results: List[Dict[str, Any]],
        output_dir: str = ".",
        all_orfs: Optional[List[Dict]] = None,
    ) -> None:
        if not isinstance(results, list):
            raise TypeError(f"results must be a list, got {type(results).__name__}")
        if not all(isinstance(r, dict) for r in results):
            raise TypeError("every item in results must be a dict")
        if not isinstance(output_dir, str) or not output_dir.strip():
            raise ValueError("output_dir must be a non-empty string")
        if all_orfs is not None:
            if not isinstance(all_orfs, list):
                raise TypeError(f"all_orfs must be a list or None, got {type(all_orfs).__name__}")
            if not all(isinstance(o, dict) for o in all_orfs):
                raise TypeError("every item in all_orfs must be a dict")

        self.results    = results
        self.output_dir = output_dir
        self.all_orfs   = all_orfs

        try:
            os.makedirs(self.output_dir, exist_ok=True)
        except OSError as e:
            raise OSError(f"Could not create output directory '{self.output_dir}': {e}") from e

    # ------------------------------------------------------------------
    # JSON
    # ------------------------------------------------------------------
    def write_json(self, filename: str = "report.json") -> None:
        """
        Purpose: Write results to a JSON file.

        Input:
           filename (str): output filename, default report.json.

        Output: JSON file.

        High-level steps:
            1. If results is empty, warn and stop.
            2. Build output path from self.output_dir + filename.
            3. Write results list to JSON using pandas.
        """
        if not self.results:
            sys.stderr.write("Warning: no results to write to JSON — skipping.\n")
            return

        if not isinstance(filename, str) or not filename.strip():
            raise ValueError("filename must be a non-empty string")

        output_path = os.path.join(self.output_dir, filename)

        try:
            df = pd.DataFrame(self.results)
            df.to_json(output_path, orient='records', indent=4)
        except (ValueError, OSError) as e:
            raise RuntimeError(f"Failed to write JSON report to '{output_path}': {e}") from e

        print(f"JSON report successfully saved to {output_path}")

    # ------------------------------------------------------------------
    # CSV
    # ------------------------------------------------------------------
    def write_csv(self, filename: str = "report.csv") -> None:
        """
        Purpose: Write results to a CSV file.

        Input:
           filename (str): output filename, default report.csv.

        Output: CSV file.

        High-level steps:
            1. If results is empty, warn and stop.
            2. Build output path from self.output_dir + filename.
            3. Extract per-sequence fields, including primary frameshift info.
            4. Write to CSV using pandas.
        """
        if not self.results:
            sys.stderr.write("Warning: no results to write to CSV — skipping.\n")
            return

        if not isinstance(filename, str) or not filename.strip():
            raise ValueError("filename must be a non-empty string")

        output_path = os.path.join(self.output_dir, filename)
        rows = []

        for i, res in enumerate(self.results):
            if not isinstance(res, dict):
                sys.stderr.write(f"Warning: result at index {i} is not a dict — skipping row.\n")
                continue

            details  = res.get("frameshift_details") or []
            primary  = details[0] if details else {}
            neighbor = primary.get("neighboring_orf") or {}

            row = {
                "sequence_id"          : res.get("sequence_id"),
                "frameshift_suspected" : res.get("frameshift_boolean"),
                "frameshift_position"  : primary.get("shift_position"),
                "shift_type"           : primary.get("shift_type"),
                "dominance"            : res.get("dominance_ratio"),
                "longest_orf_start"    : res.get("start"),
                "longest_orf_end"      : res.get("end"),
                "longest_orf_length_nt": res.get("length"),
                "overlap_length"       : neighbor.get("length"),
            }
            rows.append(row)

        if not rows:
            sys.stderr.write("Warning: no valid rows to write to CSV — skipping.\n")
            return

        cols_order = [
            "sequence_id", "frameshift_suspected", "frameshift_position",
            "shift_type", "dominance", "longest_orf_start",
            "longest_orf_end", "longest_orf_length_nt", "overlap_length",
        ]

        try:
            df   = pd.DataFrame(rows)
            cols = [c for c in cols_order if c in df.columns]
            df[cols].to_csv(output_path, index=False)
        except (ValueError, OSError) as e:
            raise RuntimeError(f"Failed to write CSV report to '{output_path}': {e}") from e

        print(f"CSV report successfully saved to {output_path}")

    # ------------------------------------------------------------------
    # HTML summary
    # ------------------------------------------------------------------
    def write_html(self, filename: str = "report.html") -> None:
        """
        Purpose: Write an HTML summary report for all sequences.

        Input:
            filename (str): output filename, default report.html.

        Output: interactive HTML dashboard.

        High-level steps:
            1. If results is empty, warn and stop.
            2. Build one bar per sequence showing backbone + longest ORF.
            3. Write HTML to output_dir.
        """
        if not self.results:
            sys.stderr.write("Warning: no results to write to HTML — skipping.\n")
            return

        if not isinstance(filename, str) or not filename.strip():
            raise ValueError("filename must be a non-empty string")

        fig = go.Figure()

        for i, res in enumerate(self.results):
            if not isinstance(res, dict):
                sys.stderr.write(f"Warning: result at index {i} is not a dict — skipping.\n")
                continue

            seq_id   = res.get("sequence_id", f"seq_{i}")
            full_len = res.get("full_seq_length", res.get("end", 0) + 100)
            length   = res.get("length")
            start    = res.get("start")
            end      = res.get("end")

            if None in (length, start, end):
                sys.stderr.write(
                    f"Warning: result '{seq_id}' missing start/end/length — skipping trace.\n"
                )
                continue

            fig.add_trace(go.Bar(
                x=[full_len], y=[seq_id], orientation='h',
                marker_color='rgba(211, 211, 211, 0.3)',
                hoverinfo='skip', showlegend=False,
            ))
            fig.add_trace(go.Bar(
                x=[length], y=[seq_id],
                base=start, orientation='h',
                name=f"Main ORF ({seq_id})",
                marker_color='mediumseagreen',
                hovertemplate=(
                    f"<b>Primary ORF</b><br>"
                    f"Coords: {start} - {end}<br>"
                    f"Length: {length}nt<extra></extra>"
                ),
            ))

        output_path = os.path.join(self.output_dir, filename)

        try:
            fig.write_html(output_path)
        except OSError as e:
            raise RuntimeError(f"Failed to write HTML report to '{output_path}': {e}") from e

        print(f"Interactive HTML dashboard successfully saved to {output_path}")

    # ------------------------------------------------------------------
    # Frameshift sequence-track plot
    # ------------------------------------------------------------------
    def write_frameshift_plot(self, filename: str = "frameshift_plot.html") -> None:
        """
        Purpose:
            Build the sequence-track plot showing the longest ORF and all
            neighboring ORF frameshift candidates across 3 reading frames,
            then write it as an interactive HTML file.

        Input:
            filename (str): output filename, default frameshift_plot.html.

        Output:
            Interactive HTML file in self.output_dir.

        High-level steps:
            1. Validate a frameshift result exists in self.results.
            2. Validate frameshift_details is present and non-empty.
            3. Validate each detail has required keys before drawing.
            4. Compute total nt from self.all_orfs if available.
            5. Draw longest ORF bar, neighbor bars, labels, vline.
            6. Zoom x-axis to 500 nt before first neighbor through last neighbor end.
            7. Write HTML to output_dir.
        """
        if not isinstance(filename, str) or not filename.strip():
            raise ValueError("filename must be a non-empty string")

        if not self.results:
            sys.stderr.write("Warning: no results provided — skipping frameshift plot.\n")
            return

        frameshift_result = next(
            (r for r in self.results if isinstance(r, dict) and r.get("frameshift_boolean")),
            None,
        )
        if not frameshift_result:
            sys.stderr.write("Warning: no result with frameshift_boolean=True — skipping plot.\n")
            return

        lo                 = frameshift_result
        frameshift_details = lo.get("frameshift_details") or []

        if not frameshift_details:
            sys.stderr.write("Warning: frameshift_details is empty — skipping plot.\n")
            return

        # Validate required keys on the longest ORF dict
        for key in ('start', 'end', 'length', 'frame', 'dominance_ratio'):
            if lo.get(key) is None:
                raise KeyError(f"Longest ORF result is missing required key '{key}'")

        if lo['frame'] not in (0, 1, 2):
            raise ValueError(f"Longest ORF frame must be 0, 1, or 2 — got {lo['frame']}")

        # Validate each detail dict before drawing
        valid_details = []
        for i, detail in enumerate(frameshift_details):
            if not isinstance(detail, dict):
                sys.stderr.write(f"Warning: frameshift_details[{i}] is not a dict — skipping.\n")
                continue
            n = detail.get("neighboring_orf")
            if not isinstance(n, dict):
                sys.stderr.write(f"Warning: frameshift_details[{i}] missing neighboring_orf — skipping.\n")
                continue
            if not all(k in n for k in ('start', 'end', 'length', 'frame')):
                sys.stderr.write(f"Warning: neighboring_orf at index {i} missing required keys — skipping.\n")
                continue
            if n['frame'] not in (0, 1, 2):
                sys.stderr.write(f"Warning: neighboring_orf at index {i} has invalid frame {n['frame']} — skipping.\n")
                continue
            if 'shift_type' not in detail:
                sys.stderr.write(f"Warning: frameshift_details[{i}] missing shift_type — skipping.\n")
                continue
            valid_details.append(detail)

        if not valid_details:
            sys.stderr.write("Warning: no valid frameshift details after validation — skipping plot.\n")
            return

        _total_nt = (
            sum(o['length'] for o in self.all_orfs if isinstance(o, dict) and 'length' in o)
            if self.all_orfs else None
        )

        LANE               = {0: 0.0, 1: 1.0, 2: 2.0}
        BAR_H              = 0.3
        COLORS             = {'longest': 'steelblue', '+1': 'darkorange', '+2': 'crimson'}
        MIN_LABEL_WIDTH_NT = 80

        fig  = go.Figure()
        y_lo = LANE[lo['frame']]

        # Longest ORF bar
        fig.add_shape(
            type='rect',
            x0=lo['start'], x1=lo['end'],
            y0=y_lo - BAR_H, y1=y_lo + BAR_H,
            fillcolor=COLORS['longest'], line_color='navy', line_width=1,
        )
        fig.add_annotation(
            x=(lo['start'] + lo['end']) / 2, y=y_lo,
            text=f"Longest ORF  |  frame {lo['frame']}  |  {lo['length']} nt  |  dom={lo['dominance_ratio']}",
            showarrow=False,
            font=dict(color='white', size=11, family='monospace'),
        )

        # Neighbor bars
        seen_labels = set()
        for detail in valid_details:
            n     = detail['neighboring_orf']
            stype = detail['shift_type']
            y     = LANE[n['frame']]
            color = COLORS.get(stype, 'grey')
            dom   = round(n['length'] / _total_nt, 3) if _total_nt else 'N/A'

            overlaps_main = (n['start'] < lo['end']) and (n['end'] > lo['start'])

            fig.add_shape(
                type='rect',
                x0=n['start'], x1=n['end'],
                y0=y - BAR_H,  y1=y + BAR_H,
                fillcolor=color, line_color='black',
                line_width=2.5 if overlaps_main else 1,
                line_dash='solid' if overlaps_main else 'dot',
                opacity=0.85 if overlaps_main else 0.55,
            )

            # Start position label above bar
            fig.add_annotation(
                x=n['start'], y=y + BAR_H + 0.05,
                xref='x', yref='y',
                text=f"nt {n['start']}",
                showarrow=False, xanchor='left',
                font=dict(color=color, size=8, family='monospace'),
            )

            # Dominance label inside bar if wide enough
            if (n['end'] - n['start']) >= MIN_LABEL_WIDTH_NT:
                fig.add_annotation(
                    x=(n['start'] + n['end']) / 2, y=y,
                    text=f"{stype}  dom={dom}",
                    showarrow=False,
                    font=dict(color='white', size=9, family='monospace'),
                )

            label = f'Neighbor ({stype}) — overlaps main' if overlaps_main else f'Neighbor ({stype})'
            if label not in seen_labels:
                seen_labels.add(label)
                fig.add_trace(go.Scatter(
                    x=[None], y=[None], mode='markers',
                    marker=dict(
                        size=12, color=color, symbol='square',
                        line=dict(color='black', width=2 if overlaps_main else 1),
                    ),
                    name=label,
                ))

        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode='markers',
            marker=dict(size=12, color=COLORS['longest'], symbol='square'),
            name=f'Longest ORF (frame {lo["frame"]}, dom={lo["dominance_ratio"]})',
        ))

        # Frameshift site vline
        shift_nt = valid_details[0]['shift_position']
        fig.add_vline(
            x=shift_nt,
            line_dash='dash', line_color='black', line_width=1.5,
            annotation_text=f'Frameshift site  nt {shift_nt}',
            annotation_position='top left',
            annotation_font=dict(color='black', size=11),
        )

        # Zoom window
        first_neighbor_start = min(d['neighboring_orf']['start'] for d in valid_details)
        last_neighbor_end    = max(d['neighboring_orf']['end']   for d in valid_details)
        x_left  = max(lo['start'], first_neighbor_start - 500)
        x_right = last_neighbor_end + 100

        fig.add_annotation(
            x=x_left, y=y_lo, xref='x', yref='y',
            text=f'← ORF start @ nt {lo["start"]}',
            showarrow=False, xanchor='left',
            font=dict(color='navy', size=10, family='monospace'),
            bgcolor='rgba(255,255,255,0.75)',
        )

        fig.update_layout(
            title='Potential Frameshift Sites — Neighboring ORFs Across 3 Reading Frames',
            xaxis=dict(title='Nucleotide Position', range=[x_left, x_right]),
            yaxis=dict(
                tickvals=list(LANE.values()),
                ticktext=[f'Frame {f}' for f in LANE],
                range=[-0.6, 2.6],
            ),
            height=450,
            plot_bgcolor='white',
            paper_bgcolor='white',
            legend=dict(x=1.01, y=1, xanchor='left'),
        )

        output_path = os.path.join(self.output_dir, filename)

        try:
            fig.write_html(output_path)
        except OSError as e:
            raise RuntimeError(f"Failed to write frameshift plot to '{output_path}': {e}") from e

        print(f"Frameshift plot successfully saved to {output_path}")
        return fig.to_html(full_html=False, include_plotlyjs='cdn')
    # ------------------------------------------------------------------
    # Dispatcher
    # ------------------------------------------------------------------
    def produce_report(self, output_format: str) -> None:
        """
        Purpose: Main entry point used by main.py to write outputs.

        Input: output_format (str): "json", "csv", "html", or "frameshift_plot".

        Output: writes a file.
        """
        if not isinstance(output_format, str) or not output_format.strip():
            raise ValueError("output_format must be a non-empty string")

        fmt = output_format.strip().lower()

        if fmt == "json":
            self.write_json()
        elif fmt == "csv":
            self.write_csv()
        elif fmt == "html":
            self.write_html()
        elif fmt == "frameshift_plot":
            self.write_frameshift_plot()
        else:
            sys.stderr.write(
                f"Unsupported report format: '{output_format}'. "
                "Use 'json', 'csv', 'html', or 'frameshift_plot'.\n"
            )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="ORF detection with frameshift analysis + sequence-track plot."
    )
    parser.add_argument("-f",  "--fasta_file",  type=str, required=True,          help="Path to input FASTA file.")
    parser.add_argument("-m",  "--min_length",  type=int, default=150,             help="Minimum ORF length (default 150).")
    parser.add_argument("-oc", "--output_csv",  type=str,                          help="Filename for CSV report.")
    parser.add_argument("-oj", "--output_json", type=str,                          help="Filename for JSON report.")
    parser.add_argument("-oh", "--output_html", type=str,                          help="Filename for HTML summary.")
    parser.add_argument("-op", "--output_plot", type=str,                          help="Filename for frameshift track plot (HTML).")
    parser.add_argument("-d",  "--output_dir",  type=str, default="../results",    help="Output directory (default ../results).")
    args = parser.parse_args()

    if not os.path.isfile(args.fasta_file):
        sys.stderr.write(f"Error: FASTA file not found: {args.fasta_file}\n")
        sys.exit(1)

    if args.min_length < 0:
        sys.stderr.write("Error: --min_length must be a positive integer.\n")
        sys.exit(1)

    if not any([args.output_csv, args.output_json, args.output_html, args.output_plot]):
        args.output_csv = "results.csv"
        sys.stdout.write("No output format specified — defaulting to CSV.\n")

    records = read_fasta(args.fasta_file)
    sys.stdout.write(f"Loaded {len(records)} sequences from {args.fasta_file}.\n")

    final_results        = []
    plot_results         = []
    all_orfs_last        = []
    all_orfs_unfiltered  = []

    for record in records:
        seq_id   = record["ID"]
        full_seq = record["Sequence"]

        # --- filtered pipeline: used for CSV / JSON / HTML ---
        all_orfs = detect_all_frames(full_seq, min_length=args.min_length)
        if not all_orfs:
            sys.stdout.write(f"  {seq_id}: no ORFs found — skipping.\n")
            continue

        detector = FrameshiftDetector(all_orfs, window=200)
        long_orf = detector.longest_orf()
        result   = detector.analyze(long_orf)
        result["sequence_id"]     = seq_id
        result["full_seq_length"] = len(full_seq)
        final_results.append(result)
        all_orfs_last = all_orfs

        # --- unfiltered pipeline: used for the plot to show all potential sites ---
        if args.output_plot:
            all_orfs_raw     = detect_all_frames(full_seq)
            det_raw          = FrameshiftDetector(all_orfs_raw, window=200)
            long_orf_raw     = det_raw.longest_orf()
            result_raw       = det_raw.analyze(long_orf_raw)
            result_raw["sequence_id"]     = seq_id
            result_raw["full_seq_length"] = len(full_seq)
            plot_results.append(result_raw)
            all_orfs_unfiltered = all_orfs_raw

        sys.stdout.write(f"  Processed: {seq_id}\n")

    if not final_results:
        sys.stdout.write("No valid ORFs found. No reports generated.\n")
        sys.exit(0)

    # reports use filtered results; plot uses unfiltered results
    report      = OrfReport(final_results,  output_dir=args.output_dir, all_orfs=all_orfs_last)
    plot_report = OrfReport(plot_results or final_results,
                            output_dir=args.output_dir,
                            all_orfs=all_orfs_unfiltered or all_orfs_last)

    if args.output_csv:
        report.write_csv(args.output_csv)
    if args.output_json:
        report.write_json(args.output_json)
    if args.output_html:
        report.write_html(args.output_html)
    if args.output_plot:
        plot_report.write_frameshift_plot(args.output_plot)
