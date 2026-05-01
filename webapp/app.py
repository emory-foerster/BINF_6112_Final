#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from flask import Flask, send_from_directory, abort

APP_DIR = os.path.dirname(os.path.abspath(__file__))
CSS_DIR = os.path.join(APP_DIR, "css")
ROOT_DIR = os.path.dirname(APP_DIR)
OUT_DIR = os.path.join(ROOT_DIR, "out_results")
REPORT_HTML = "full_report.html"
FRAMESHIFT_HTML = "frameshift_plot.html"

app = Flask(__name__, static_folder=None)

@app.get("/")
def home():
    return send_from_directory(APP_DIR, "index.html")

@app.get("/css/<path:filename>")
def css(filename):
    return send_from_directory(CSS_DIR, filename)

@app.get("/report")
def report():
    path = os.path.join(OUT_DIR, REPORT_HTML)
    if not os.path.isfile(path):
        abort(404, description=f"Missing {REPORT_HTML} in out_results/")

    return send_from_directory(OUT_DIR, REPORT_HTML)

@app.get("/frameshift-plot")
def frameshift_plot():
    path = os.path.join(OUT_DIR, FRAMESHIFT_HTML)
    if not os.path.isfile(path):
        abort(404, description=f"Missing {FRAMESHIFT_HTML} in out_results/")
    return send_from_directory(OUT_DIR, FRAMESHIFT_HTML)

@app.get("/files")

def files():
    return "<h1>404. This page is unavailable right now</h1>"

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5500, debug=True)