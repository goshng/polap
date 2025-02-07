import os
import subprocess
import pandas as pd


def run_command(command):
    """Run a shell command and handle errors."""
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
        exit(1)


def blast_alignment(seq1, seq2, output_file, outfmt, task="blastn", strand="both"):
    """Run BLAST to align two sequences and save results."""
    command = f'blastn -task {task} -strand {strand} -query {seq1} -subject {seq2} -out {output_file} -outfmt "{outfmt}"'
    run_command(command)


def parse_blast_output(blast_output, output_dir, seq2_length):
    """Parse the BLAST output and return the relevant qstart or adjusted sstart value."""
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


def get_sequence_length(fasta_file):
    """Get the length of the sequence in the FASTA file."""
    result = subprocess.run(
        f"seqkit stats {fasta_file} --tabular | tail -n 1",
        shell=True,
        capture_output=True,
        text=True,
    )
    return int(result.stdout.split()[4])  # Extract length from seqkit output


def compute_coverage1(blast_output, seq1_length, output_dir):
    """Compute the alignment coverage based on the merged bed."""
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


def compute_coverage2(blast_output, seq2_length, output_dir):
    """Compute the alignment coverage based on the merged bed."""
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


def get_pident(blast_output, output_dir):
    """
    Read the BLAST output and check the first line for specific conditions.
    Write the pident value to a file in the output directory. If conditions are not met, write 0.
    """
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


def main(seq1, seq2, output_dir):
    """Main workflow to determine alignment of two sequences by restart."""
    os.makedirs(output_dir, exist_ok=True)

    # Get the length of seq1
    seq1_length = get_sequence_length(seq1)

    # Get the length of seq2
    seq2_length = get_sequence_length(seq2)

    # First BLAST alignment
    blast_output_1 = os.path.join(output_dir, "blast_results_1.txt")
    print(f"Running first BLAST alignment between {seq1} and {seq2}...")
    blast_alignment(
        seq1,
        seq2,
        blast_output_1,
        "6 qseqid qstart qend sseqid sstart send sstrand pident length qlen slen nident mismatch gaps",
    )

    # Parse the BLAST output to find the relevant position
    relevant_position = parse_blast_output(blast_output_1, output_dir, seq2_length)

    print(f"Using relevant position {relevant_position} to restart {seq2}.")

    # Restart seq2 using the relevant position
    seq2_restarted = os.path.join(output_dir, "seq2_restarted.fasta")
    run_command(f"seqkit restart -i {relevant_position} {seq2} -o {seq2_restarted}")

    # Second BLAST alignment with dc-megablast and strand plus
    blast_output_2 = os.path.join(output_dir, "blast_results_2.txt")
    print(
        f"Running second BLAST alignment with restarted {seq2_restarted} using dc-megablast and strand plus..."
    )
    blast_alignment(
        seq1,
        seq2_restarted,
        blast_output_2,
        "6 qseqid qstart qend sseqid sstart send sstrand pident length qlen slen nident mismatch gaps",
        task="dc-megablast",
        strand="plus",
    )

    # Compute the alignment coverage
    coverage1 = compute_coverage1(blast_output_2, seq1_length, output_dir)
    coverage2 = compute_coverage2(blast_output_2, seq2_length, output_dir)
    print(f"Alignment coverage (seq1-based): {coverage1}")
    print(f"Alignment coverage (seq2-based): {coverage2}")

    print(f"Second BLAST results saved to {blast_output_2}")

    try:
        pident = get_pident(blast_output_2, output_dir)
        print(f"The pident value is saved to {output_dir}/pident.txt and is: {pident}")
    except ValueError as e:
        print(e)


if __name__ == "__main__":
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

    main(args.seq1, args.seq2, args.out)
