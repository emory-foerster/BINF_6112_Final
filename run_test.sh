#!/bin/bash
#SBATCH --job-name=frame_report
#SBATCH --partition=Centaurus
#SBATCH --cpus-per-task=8
#SBATCH --time=10:30:00
#SBATCH --mem-per-cpu=2G
#SBATCH --output=job_%j.log

# run_test.sh
# Validates the full ORF-frameshift pipeline end-to-end.
# Run from the repo root: bash run_test.sh
#
# Flags used:
#   -f   Input FASTA (query sequences)
#   -r   Reference FASTA (single sequence, used for comparison baseline)
#   -m   Minimum ORF length in bp (200 filters noise)
#   -d   Output directory
#   -oc  CSV report filename
#   -oj  JSON report filename
#   -oh  Full interactive HTML report (html_report2.py)

# Resolve the repo root from wherever this script lives
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INPUT_FILE="$REPO_ROOT/datasets/seq1.fasta"
REF_FILE="$REPO_ROOT/datasets/sequence_ref.fasta"
OUTPUT_DIR="$REPO_ROOT/out_results"
CSV_OUT="test_summary.csv"
JSON_OUT="test_data.json"
HTML_OUT="test_report.html"

mkdir -p "$OUTPUT_DIR"

echo "Starting ORF-Frameshift Pipeline Test..."
echo "Output directory: $OUTPUT_DIR"
echo "Running main.py ..."

python "$REPO_ROOT/src/main.py" \
  -f "$INPUT_FILE" \
  -r "$REF_FILE" \
  -m 200 \
  -d "$OUTPUT_DIR" \
  -oc "$CSV_OUT" \
  -oj "$JSON_OUT" \
  -oh "$HTML_OUT" \
  > "$OUTPUT_DIR/frame_job_out.txt" \
  2> "$OUTPUT_DIR/frame_job_err.txt"

echo "Pipeline finished. Checking outputs..."

for OUTFILE in "$OUTPUT_DIR/$CSV_OUT" "$OUTPUT_DIR/$JSON_OUT" "$OUTPUT_DIR/$HTML_OUT"; do
  if [ -f "$OUTFILE" ]; then
    echo "  SUCCESS: $OUTFILE"
  else
    echo "  ERROR:   $OUTFILE missing"
  fi
done

echo ""
echo "Stdout log: $OUTPUT_DIR/frame_job_out.txt"
echo "Stderr log: $OUTPUT_DIR/frame_job_err.txt"
echo "Test Complete."
