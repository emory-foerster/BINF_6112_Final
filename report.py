"""
report.py

Writes the final results to an output file.
Supports CSV or JSON output, based on the user's chosen format.
Used as the last step of the pipeline called by main.py.
"""

def produce_report(format):
    """
    1. If format is "csv":
            write a header row using the keys
            write one row per result
    2. If format is "json":
            write the list of dicts as JSON
    3. If format is not supported:
            return an error message.
    """
    pass