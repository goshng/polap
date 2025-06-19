################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################

# polaplib/run-polap-py-compare2ptdna.py

################################################################################
# We use BLAST to roughly align two very similar sequences if they are
# only different in their start position. 
# Of the 4 possible ptDNA sequences, one of them must be aligned to a
# reference sequence because they are very similar. 
# So, we roughly try to find the new start position of the second sequence,
# so that we can globally align the first sequence.
# To do this, we BLAST twice; first to find the start position and
# second to check if the first and the restarted second sequences are
# well aligned.
# To check this 2nd part, we check only the top query and subject alignment
# in the BLAST. If the start and end of the query and those of the subject
# are within 1% of each sequence length, then we use the pident value
# to report. This works fine because we use four such values from the
# 4 possible ptDNA in order to choose one of them as the one similar to
# the reference.
# Align two ptDNA sequences using seqkit restart and BLAST.
#
# Example:
# python polap-py-compare2ptdna.py \
#   --seq1 input/ptdna-reference.fa \
#   --seq2 input/run-polap-py-compare2ptdna.input1.5.fa \
#   --out output/run-polap-py-compare2ptdna.output1.5.52-mtdna/1
#
# This script compares two ptDNA sequences.
#
# Used by:
# function _run_polap_compare2ptdna {
# function _disassemble-step14 {
# function _disassemble-step15 {
#
################################################################################

# %%
import os
import subprocess
import pandas as pd
import argparse

debug = os.getenv("_POLAP_DEBUG", "0")
# print(f"[script.py] DEBUG = {debug}")

# %%
is_test = True
is_test = False

# %%
def run_command(command):
    """Run a shell command and handle errors."""
    """Comment"""
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
        exit(1)


# %%
# BLAST on two sequences.
def blast_alignment(seq1, seq2, output_file, outfmt, task="blastn", strand="both"):
    command = f'blastn -task {task} -strand {strand} -query {seq1} -subject {seq2} -out {output_file} -outfmt "{outfmt}"'
    run_command(command)


# %%
# Parse the BLAST output and return the relevant qstart or adjusted sstart value.
def parse_blast_output(blast_output, output_dir, seq2_length):
    columns = [
        "qseqid",
        "qstart",
        "qend",
        "sseqid",
        "sstart",
        "send",
        "sstrand",
        "pident",
        "length",
        "qlen",
        "slen",
        "nident",
        "mismatch",
        "gaps",
    ]
    df = pd.read_csv(blast_output, sep="\t", names=columns)

    # Case 1: Locate the row where qstart == 1 and sstrand == "plus"
    for i in range(1, 10):
        match_qstart = df[(df["qstart"] == i) & (df["sstrand"] == "plus")]
        if not match_qstart.empty:
            sstart = match_qstart.iloc[0]["sstart"]
            qstart = match_qstart.iloc[0]["qstart"]
            if qstart <= sstart:
                adjusted_value = sstart - (qstart - 1)
                return adjusted_value
            # return match_qstart.iloc[0]["sstart"]

    # Case 2: Locate the row where sstart == 1 and sstrand == "plus"
    for i in range(1, 10):
        match_qstart = df[(df["sstart"] == i) & (df["sstrand"] == "plus")]
        if not match_qstart.empty:
            sstart = match_qstart.iloc[0]["sstart"]
            qstart = match_qstart.iloc[0]["qstart"]
            if qstart > sstart:
                adjusted_value = seq2_length - (qstart - i)
                return adjusted_value

    # Case 3: Neither case found, handle error
    print("Error: No valid match found in BLAST output.")
    coverage_file = os.path.join(output_dir, "coverage.txt")
    with open(coverage_file, "w") as f:
        f.write("0\n")
    exit(1)


# %%
# Get the length of the sequence in the FASTA file.
def get_sequence_length(fasta_file):
    result = subprocess.run(
        f"seqkit stats {fasta_file} --tabular | tail -n 1",
        shell=True,
        capture_output=True,
        text=True,
    )
    # Extract length from seqkit output
    return int(result.stdout.split()[4])


# %%
# Compute the alignment coverage based on the merged bed.
def compute_coverage1(blast_output, seq1_length, output_dir):
    bed_file = os.path.join(output_dir, "blast_results_2.bed")
    merged_bed_file = os.path.join(output_dir, "merged1.bed")

    # Extract the first 5 lines and convert to BED format
    run_command(
        f"head -n 25 {blast_output} | cut -f1-3 | sort -k1,1 -k2,2n > {bed_file}"
    )

    # Merge the BED file and compute total covered length
    coverage_output = subprocess.run(
        f"bedtools merge -i {bed_file} | awk -F'\t' 'BEGIN{{SUM=0}}{{ SUM+=$3-$2 }}END{{print SUM}}'",
        shell=True,
        capture_output=True,
        text=True,
    )
    if coverage_output.stdout.strip():
        covered_length = int(coverage_output.stdout.strip())
    else:
        covered_length = 0

    # Compute coverage
    coverage1 = covered_length / seq1_length

    # Save coverage to file
    coverage_file = os.path.join(output_dir, "coverage1.txt")
    with open(coverage_file, "w") as f:
        f.write(f"{coverage1}\n")

    return coverage1


# %%
# Compute the alignment coverage based on the merged bed.
def compute_coverage2(blast_output, seq2_length, output_dir):
    bed_file = os.path.join(output_dir, "blast_results_2.bed")
    merged_bed_file = os.path.join(output_dir, "merged2.bed")

    # Extract the first 5 lines and convert to BED format
    run_command(
        f"head -n 25 {blast_output} | cut -f4-6 | sort -k1,1 -k2,2n > {bed_file}"
    )

    # Merge the BED file and compute total covered length
    coverage_output = subprocess.run(
        f"bedtools merge -i {bed_file} | awk -F'\t' 'BEGIN{{SUM=0}}{{ SUM+=$3-$2 }}END{{print SUM}}'",
        shell=True,
        capture_output=True,
        text=True,
    )
    if coverage_output.stdout.strip():
        covered_length = int(coverage_output.stdout.strip())
    else:
        covered_length = 0

    # Compute coverage
    coverage2 = covered_length / seq2_length

    # Save coverage to file
    coverage_file = os.path.join(output_dir, "coverage2.txt")
    with open(coverage_file, "w") as f:
        f.write(f"{coverage2}\n")

    return coverage2


# %%
# Read the BLAST output and check the first line for specific conditions.
# Write the pident value to a file in the output directory. If conditions are not met, write 0.
def get_pident(blast_output, output_dir):
    # Define column names for the BLAST output
    columns = [
        "qseqid",
        "qstart",
        "qend",
        "sseqid",
        "sstart",
        "send",
        "sstrand",
        "pident",
        "length",
        "qlen",
        "slen",
        "nident",
        "mismatch",
        "gaps",
    ]

    # Output file path
    pident_file = os.path.join(output_dir, "pident.txt")

    # Read the BLAST output into a DataFrame
    try:
        df = pd.read_csv(blast_output, sep="\t", names=columns)
    except Exception as e:
        raise ValueError(f"Error reading the BLAST file: {e}")

    # Ensure the file has data
    if df.empty:
        with open(pident_file, "w") as f:
            f.write("0\n")
        return 0

    # Check the conditions on the first line
    first_row = df.iloc[0]
    if (
        first_row["qstart"] < first_row["qlen"] * 0.01
        and first_row["qend"] > first_row["qlen"] * 0.99
        and first_row["sstart"] < first_row["slen"] * 0.01
        and first_row["send"] > first_row["slen"] * 0.99
    ):
        pident = first_row["pident"]
    else:
        pident = 0

    # Write the pident value to the output file
    os.makedirs(output_dir, exist_ok=True)
    with open(pident_file, "w") as f:
        f.write(f"{pident}\n")

    return pident


# %%
# Main workflow to determine alignment of two sequences by restart.
def main(seq1, seq2, output_dir):

# %%
    if is_test:
        seq1 = args.seq1
        seq2 = args.seq2
        output_dir = args.out
        print("Test is on.")

# %%
    os.makedirs(output_dir, exist_ok=True)

# %%
    # Get the length of seq1
    seq1_length = get_sequence_length(seq1)

    # Get the length of seq2
    seq2_length = get_sequence_length(seq2)

# %%
    # First BLAST alignment
    blast_output_1 = os.path.join(output_dir, "blast_results_1.txt")
    print(f"Running first BLAST alignment between {seq1} and {seq2}...")
    blast_alignment(
        seq1,
        seq2,
        blast_output_1,
        "6 qseqid qstart qend sseqid sstart send sstrand pident length qlen slen nident mismatch gaps",
    )

# %%
    # Parse the BLAST output to find the relevant position
    relevant_position = parse_blast_output(blast_output_1, output_dir, seq2_length)
    print(f"Using relevant position {relevant_position} to restart {seq2}.")

# %%
    # Restart seq2 using the relevant position
    seq2_restarted = os.path.join(output_dir, "seq2_restarted.fasta")
    run_command(f"seqkit restart -i {relevant_position} {seq2} -o {seq2_restarted}")

# %%
    # Second BLAST alignment with dc-megablast and strand plus
    blast_output_2 = os.path.join(output_dir, "blast_results_2.txt")
    print(
        f"Running second BLAST alignment with restarted {seq2_restarted}"
        f" using dc-megablast and strand plus..."
    )
    blast_alignment(
        seq1,
        seq2_restarted,
        blast_output_2,
        "6 qseqid qstart qend sseqid sstart send sstrand pident length qlen slen nident mismatch gaps",
        task="dc-megablast",
        strand="plus",
    )

# %%
    # Compute the alignment coverage
    coverage1 = compute_coverage1(blast_output_2, seq1_length, output_dir)
    coverage2 = compute_coverage2(blast_output_2, seq2_length, output_dir)
    print(f"Alignment coverage (seq1-based): {coverage1}")
    print(f"Alignment coverage (seq2-based): {coverage2}")

    print(f"Second BLAST results saved to {blast_output_2}")

# %%
    try:
        pident = get_pident(blast_output_2, output_dir)
        print(f"The pident value is saved to {output_dir}/pident.txt and is: {pident}")
    except ValueError as e:
        print(e)

# %%
################################################################################
# Tests
def test_f1():

# %%
    parser = argparse.ArgumentParser(
        description="Align two ptDNA sequences using seqkit restart and BLAST."
    )
    parser.add_argument(
        "--seq1",
        required=True,
        help="Path to the first ptDNA sequence in FASTA format.",
    )
    parser.add_argument(
        "--seq2",
        required=True,
        help="Path to the second ptDNA sequence in FASTA format.",
    )
    parser.add_argument(
        "--out", required=True, help="Output directory to store results."
    )

    # Test arguments
    custom_args = [
        '--seq1', 'input/ptdna-reference.fa', 
        '--seq2', 'input/run-polap-py-compare2ptdna.input1.5.fa', 
        '--out', 'output/run-polap-py-compare2ptdna.output1.5.52-mtdna/1'
    ]
    args = parser.parse_args(custom_args)

# %%
    print(args)

# %%
    seq1 = args.seq1
    seq2 = args.seq2
    output_dir = args.out
    print("Test is on.")
    print(f"seq1: {seq1}")
    print(f"seq2: {seq2}")
    print(f"output_dir: {output_dir}")

# %%
    os.makedirs(output_dir, exist_ok=True)
    seq1_length = get_sequence_length(seq1)
    seq2_length = get_sequence_length(seq2)
    print(f"seq1 len: {seq1_length}")
    print(f"seq2 len: {seq2_length}")

# %%
    blast_output_1 = os.path.join(output_dir, "blast_results_1.txt")
    print(f"blast output1: {blast_output_1}")
    blast_alignment(
        seq1,
        seq2,
        blast_output_1,
        "6 qseqid qstart qend sseqid sstart send sstrand pident length qlen slen nident mismatch gaps",
    )

# %%
    blast_output = blast_output_1
    columns = [
        "qseqid",
        "qstart",
        "qend",
        "sseqid",
        "sstart",
        "send",
        "sstrand",
        "pident",
        "length",
        "qlen",
        "slen",
        "nident",
        "mismatch",
        "gaps",
    ]
    df = pd.read_csv(blast_output, sep="\t", names=columns)

    # Because two sequences are similar, we expect the query front
    # part would be aligned to either some of the front or back of
    # the subject sequence. Using this difference, we adjust the
    # 2nd subject sequence to restart in a different position.
    # This would make the two sequences could be globally better
    # aligned in a 2nd round of BLAST.
    # We consider only the plus alignment because ptDNA has 4
    # possible sequences to compare and one of them is expected to
    # be aligned globally. We ignore the minus.
    #
    # Case 1: Locate the row where qstart == 1 and sstrand == "plus"
    for i in range(1, 10):
        match_qstart = df[(df["qstart"] == i) & (df["sstrand"] == "plus")]
        if not match_qstart.empty:
            sstart = match_qstart.iloc[0]["sstart"]
            qstart = match_qstart.iloc[0]["qstart"]
            if qstart <= sstart:
                adjusted_value = sstart - (qstart - 1)
                # return adjusted_value
            # return match_qstart.iloc[0]["sstart"]

    # Case 2: Locate the row where sstart == 1 and sstrand == "plus"
    for i in range(1, 10):
        match_qstart = df[(df["sstart"] == i) & (df["sstrand"] == "plus")]
        if not match_qstart.empty:
            sstart = match_qstart.iloc[0]["sstart"]
            qstart = match_qstart.iloc[0]["qstart"]
            if qstart > sstart:
                adjusted_value = seq2_length - (qstart - i)
                # return adjusted_value
# %%
    relevant_position = adjusted_value
    print(f"restart position: {relevant_position}")

# %%
    # Restart seq2 using the relevant position
    seq2_restarted = os.path.join(output_dir, "seq2_restarted.fasta")
    run_command(f"seqkit restart -i {relevant_position} {seq2} -o {seq2_restarted}")

# %%
    # Second BLAST alignment with dc-megablast and strand plus
    blast_output_2 = os.path.join(output_dir, "blast_results_2.txt")
    print(
        f"Running second BLAST alignment with restarted {seq2_restarted}"
        f" using dc-megablast and strand plus..."
    )
    blast_alignment(
        seq1,
        seq2_restarted,
        blast_output_2,
        "6 qseqid qstart qend sseqid sstart send sstrand pident length qlen slen nident mismatch gaps",
        task="dc-megablast",
        strand="plus",
    )

# %%
    # Compute the alignment coverage
    coverage1 = compute_coverage1(blast_output_2, seq1_length, output_dir)
    coverage2 = compute_coverage2(blast_output_2, seq2_length, output_dir)
    print(f"Alignment coverage (seq1-based): {coverage1}")
    print(f"Alignment coverage (seq2-based): {coverage2}")

    print(f"Second BLAST results saved to {blast_output_2}")

# %%
    pident = get_pident(blast_output_2, output_dir)
    print(f"The pident value is saved to {output_dir}/pident.txt and is: {pident}")

################################################################################
# Main
# %%

if __name__ == "__main__":

# %%
    import argparse

    parser = argparse.ArgumentParser(
        description="Align two ptDNA sequences using seqkit restart and BLAST."
    )
    parser.add_argument(
        "--seq1",
        required=True,
        help="Path to the first ptDNA sequence in FASTA format.",
    )
    parser.add_argument(
        "--seq2",
        required=True,
        help="Path to the second ptDNA sequence in FASTA format.",
    )
    parser.add_argument(
        "--out", required=True, help="Output directory to store results."
    )

    args = parser.parse_args()

# %%
    main(args.seq1, args.seq2, args.out)
