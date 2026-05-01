import pytest
import os
import json
import sys
from fasta_io import read_fasta
from orf import detect_all_frames
from report import OrfReport
from main import main as run_main

# TEST 1: FASTA Parsing (Positive)
def test_fasta_parsing(tmp_path):
    """Ensures read_fasta parses multiple records correctly."""
    fasta_content = ">Seq1\nATGCGTGTC\n>Seq2 Severe acute respiratory syndrome coronavirus 2\nTAGCTATACG"
    d = tmp_path / "data"
    d.mkdir()
    f = d / "test.fasta"
    f.write_text(fasta_content)

    records = read_fasta(str(f))
    
    assert len(records) == 2
    assert records[0]["ID"] == "Seq1"
    assert records[1]["Description"] == "Severe acute respiratory syndrome coronavirus 2"

# TEST 2: ORF Detection (Positive)
def test_orf_detection():
    """Confirms an ORF is found in a simple sequence."""
    test_seq = "ATGGAGAGCCTTTAA" # Length 15
    # detect_all_frames returns a list of dictionaries
    orfs = detect_all_frames(test_seq, min_length=0)
    
    assert len(orfs) >= 1
    # Check the attributes of the first ORF found
    assert orfs[0]["start"] == 0
    assert orfs[0]["length"] == 15

# TEST 3: JSON Report Integrity (Positive)
def test_json_report_writer(tmp_path):
    """Confirms reporting logic creates valid, loadable JSON."""
    results = [{
        "sequence_id": "TestSeq",
        "frameshift_boolean": True,
        "dominance_ratio": 0.5,
        "start": 10,
        "end": 20,
        "length": 10
    }]
    output_dir = str(tmp_path / "results")
    filename = "test_report.json"
    
    # Initialize the report engine
    report = OrfReport(results, output_dir=output_dir)
    report.write_json(filename)
    
    file_path = os.path.join(output_dir, filename)
    assert os.path.exists(file_path)
    
    with open(file_path, "r") as f:
        data = json.load(f)
    
    assert data[0]["sequence_id"] == "TestSeq"