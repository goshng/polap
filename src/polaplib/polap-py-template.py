#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Help me write a template python script that I copy to kick-starting python
# coding. Prepend import statements. I would use file operation, command-line
# processing, processing execution, data processing, and the like.
# This may be a useful for bioinformatics, logging, data science.
#
# Requirements:
# biopython pandas matplotlib seaborn

# ---------- Standard Libraries ----------
import os
import sys
import argparse
import logging
import subprocess
import shutil
import json
import csv
import re
from pathlib import Path
from typing import List, Dict, Optional

# ---------- Bioinformatics Libraries ----------
from Bio import SeqIO  # Biopython for FASTA/FASTQ parsing

# ---------- Data Science Libraries ----------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# ---------- Logging Configuration ----------
def setup_logger(name: str = "bio_logger", level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(level)
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)s - %(message)s", "%H:%M:%S"
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


logger = setup_logger()

# ---------- Helper Functions ----------


def run_command(
    cmd: List[str], capture_output: bool = False, check: bool = True
) -> subprocess.CompletedProcess:
    logger.info(f"Running command: {' '.join(cmd)}")
    return subprocess.run(cmd, capture_output=capture_output, check=check, text=True)


def read_fasta(filepath: str) -> List[SeqIO.SeqRecord]:
    """Read FASTA sequences using Biopython."""
    logger.debug(f"Reading FASTA from {filepath}")
    return list(SeqIO.parse(filepath, "fasta"))


def read_fastq(filepath: str) -> List[SeqIO.SeqRecord]:
    """Read FASTQ sequences using Biopython."""
    logger.debug(f"Reading FASTQ from {filepath}")
    return list(SeqIO.parse(filepath, "fastq"))


def save_dataframe(df: pd.DataFrame, path: str) -> None:
    logger.info(f"Saving DataFrame to {path}")
    df.to_csv(path, sep="\t", index=False)


def ensure_dir(path: str) -> None:
    logger.debug(f"Ensuring directory exists: {path}")
    Path(path).mkdir(parents=True, exist_ok=True)


# ---------- Argument Parsing ----------


def parse_args():
    parser = argparse.ArgumentParser(
        description="Template for bioinformatics + data science Python scripts."
    )
    parser.add_argument("input", help="Input file or directory (FASTA/FASTQ/TSV/CSV)")
    parser.add_argument(
        "-o", "--output", help="Output file/directory", default="output/"
    )
    parser.add_argument(
        "-l", "--log", help="Log level (DEBUG, INFO, WARNING)", default="INFO"
    )
    return parser.parse_args()


# ---------- Main Logic ----------


def main():
    args = parse_args()

    # Set log level
    logger.setLevel(getattr(logging, args.log.upper(), logging.INFO))

    logger.info("Script started.")
    logger.debug(f"Args: {args}")

    ensure_dir(args.output)

    ext = Path(args.input).suffix.lower()

    if ext in [".fasta", ".fa"]:
        records = read_fasta(args.input)
        logger.info(f"Read {len(records)} FASTA records.")
    elif ext in [".fastq", ".fq"]:
        records = read_fastq(args.input)
        logger.info(f"Read {len(records)} FASTQ records.")
    elif ext in [".csv", ".tsv"]:
        sep = "\t" if ext == ".tsv" else ","
        df = pd.read_csv(args.input, sep=sep)
        logger.info(f"Read dataframe with shape {df.shape}.")
        # Example data processing
        if "length" in df.columns:
            logger.debug("Plotting length distribution.")
            sns.histplot(df["length"], kde=True)
            plt.title("Length Distribution")
            plt.xlabel("Length")
            plt.ylabel("Frequency")
            plt.savefig(Path(args.output) / "length_distribution.png")
            plt.close()
            logger.info("Plot saved.")
    else:
        logger.error(f"Unsupported file extension: {ext}")
        sys.exit(1)

    logger.info("Script finished.")


# ---------- Entry Point ----------
if __name__ == "__main__":
    main()
