#!/bin/bash

# run_test.sh - Automated Pipeline Validation
echo "Starting ORF-Frameshift Pipeline Test..."


# 1. Define paths
INPUT_FILE="../datasets/Covid_GCF_009858895.2/sample.fasta"
OUTPUT_DIR="../out_results"
CSV_OUT="test_summary.csv"
JSON_OUT="test_data.json"

# 2. Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# 3. Execute the pipeline
echo "Running main.py ..."
python main.py -f $INPUT_FILE -oc $CSV_OUT -oj $JSON_OUT -d $OUTPUT_DIR -m 200

# 4. Check if files were created
if [ -f "$OUTPUT_DIR/$CSV_OUT" ]; then
    echo "SUCCESS CSV Report generated."
else
    echo "ERROR CSV Report missing."
fi

if [ -f "$OUTPUT_DIR/$JSON_OUT" ]; then
    echo "SUCCESS JSON Data generated."
else
    echo "ERROR JSON Data missing."
fi

echo "Test Sequence Complete."