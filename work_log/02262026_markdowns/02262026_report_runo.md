# report.py (Work Summary)

## What was implemented
This file is responsible for writing the final project results to disk. For now, the report functionality supports **JSON output only**. The script takes a list of per-sequence result dictionaries and saves them into a JSON file using `json.dump()` with indentation for readability.

## Functions implemented

### `write_json(results, output_path)`

**Input:**  
- `results`: a list of dictionaries, where each dictionary contains the analysis results for one sequence  
- `output_path`: the filename to write (example: `"report.json"`)

**Output:**  
- Creates a JSON file at `output_path` containing the list of result dictionaries

### `produce_report(results, output_path, output_format)`
**Purpose:**  
Main entry point for report generation. Chooses the correct writer based on the requested format.

**Input:**  
- `results`: list of result dictionaries  
- `output_path`: output filename  
- `output_format`: format string (currently only `"json"` is supported)

**Output:**  
- Calls `write_json()` and writes the JSON report file  
## Example input
results = [
    {
        "sequence_id": "seqA",
        "frameshift_suspected": True,
        "frameshift_position": 600,
        "shift_type": "+1",
        "dominance": 0.78,
        "longest_orf_start": 120,
        "longest_orf_end": 1320,
        "longest_orf_length_nt": 1200
    }
]