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

# polaplib/run-polap-py-find-plastid-gfa2fasta.py

################################################################################
# We extract DNA sequences of simple circles with at least one edge that
# belongs to the given seeds.
# Searching the graph for all circular paths can take too long mostly
# because of no organelle genome assembly case. We use number 30 to
# stop the search.
# 
# See Also: polap-py-find-cc-with-seeds.py
#
# Example:
# mkdir output/run-polap-py-find-plastid-gfa2fasta.output.5.52-mtdna
# python polap-py-find-plastid-gfa2fasta.py \
#   --gfa input/run-polap-py-find-plastid-gfa2fasta.input1.5.gfa \
#   --seed input/run-polap-py-find-plastid-gfa2fasta.input2.5.mt.contig.name \
#   --out output/run-polap-py-find-plastid-gfa2fasta.output.5.52-mtdna
#
# input/run-polap-py-find-plastid-gfa2fasta.input1.5.gfa
# ------------------------------------------------------
# L       edge_1  +       edge_4  +       0M      RC:i:59
# L       edge_1  -       edge_4  +       0M      RC:i:55
# L       edge_2  +       edge_4  -       0M      RC:i:69
# L       edge_2  -       edge_4  -       0M      RC:i:46
# P       contig_3        edge_3+ *
# P       contig_1        edge_1+ *
# P       contig_2        edge_2+ *
# P       contig_4        edge_4+ *
#
# input/run-polap-py-find-plastid-gfa2fasta.input2.5.mt.contig.name
# -----------------------------------------------------------------
# edge_1
# edge_2
# edge_4
#
# output/run-polap-py-find-plastid-gfa2fasta.output.5.52-mtdna
# ------------------------------------------------------------
# circular_path_1_concatenated.fa
# circular_path_1.fa
# circular_path_2_concatenated.fa
# circular_path_2.fa
# circular_path_3_concatenated.fa
# circular_path_3.fa
# circular_path_4_concatenated.fa
# circular_path_4.fa
# circular_path_count.txt
# circular_path_nodes.txt
# circular_path.txt
#
# This script is used to parse gfa to find plastid DNA sequences.
#
# Used by:
# function _run_polap_disassemble {
# function _disassemble-step13 {
# function _run_polap_ptgaul {
#
################################################################################

# %%
import os
import subprocess
import sys
import networkx as nx
import argparse

debug = os.getenv("_POLAP_DEBUG", "0")

# %%
is_test = True
is_test = False

# %%
# We add two edges from one using the reverse complement.
# We need this to find all sequences from a ptDNA-graph.
# Parse a GFA file and construct a directed graph, 
# generating reverse complement edges.
def parse_gfa_to_graph(gfa_file):

# %%
    if is_test:
        gfa_file = args.gfa
        print("Test is on.")

# %%
    G = nx.DiGraph()

# %%
    with open(gfa_file, "r") as file:
        for line in file:
            if line.startswith("L"):  # Process Link lines
                fields = line.strip().split("\t")
                if len(fields) >= 5:
                    source, source_orient, target, target_orient = fields[1:5]

                    # Add the original edge
                    G.add_edge(f"{source}{source_orient}", f"{target}{target_orient}")

                    # Add the reverse complement edge
                    rev_source_orient = "-" if source_orient == "+" else "+"
                    rev_target_orient = "-" if target_orient == "+" else "+"
                    G.add_edge(
                        f"{target}{rev_target_orient}", f"{source}{rev_source_orient}"
                    )

# %%
    return G


# %%
# Find all circular paths, including single-node cycles (self-loops).
def find_circular_paths(G):
    simple_cycles = list(nx.simple_cycles(G))
    return simple_cycles


# %%
# Filter circular paths to keep only those that contain at least one node 
# from the seed file. 
# The seed contanis edge sequences with organelle genes.
#
# Use:
# filtered_circular_paths = filter_circular_paths_by_seed(
#     all_circular_paths, args.seed
# )
def filter_circular_paths_by_seed(circular_paths, seed_file):

# %%
    if is_test:
        circular_paths = all_circular_paths
        seed_file = args.seed
        print("Test is on.")

# %%
    if not seed_file:
        return circular_paths

# %%
    # Load seed nodes into a set
    with open(seed_file, "r") as f:
        seed_nodes = {line.strip() for line in f if line.strip()}

# %%
    # Filter circular paths
    # circular_paths = [
    #     ["A+", "B-", "C+"],
    #     ["X-", "Y+", "Z-"],
    #     ["M+", "N-", "O+"]
    # ]
    # 
    # seed_nodes = {"A", "Y", "Z"}
    # filtered_paths = [
    #     ["A+", "B-", "C+"],     # contains "A"
    #     ["X-", "Y+", "Z-"]      # contains "Y" and "Z"
    # ]
    filtered_paths = [
        path
        for path in circular_paths
        if any(node.rstrip("+-") in seed_nodes for node in path)
    ]

# %%
    return filtered_paths


# %%
# Run a shell command.
# Exit if we have errors.
def run_command(command):
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
        exit(1)


# %%
# Extract a sequence for a given edge with the correct orientation.
#
# Use:
# sequence = extract_edge_sequence(fasta_file, edge, "+")
def extract_edge_sequence(fasta_file, edge, orientation):

# %%
    if is_test:
        # fasta_file = given the same
        edge = "edge_4"
        orientation = "+"
        print("Test is on.")

# %%
    # A better way
    # base_edge = edge.rstrip("+-")
    base_edge = edge.split("+")[0].split("-")[0]  # Remove signs to get base edge name

# %%
    # These two are the same; one or multiline strings notation.
    # command = (
    #     f"line 1 "
    #     f"line 2"
    # )
    #
    # command = "line 1 line 2"
    if orientation == "+":
        command = f"seqkit grep -p {base_edge} {fasta_file} | seqkit seq -t dna -v"
    elif orientation == "-":
        command = (
            f"seqkit grep -p {base_edge} {fasta_file} | seqkit seq -t dna -v -r -p"
        )
    else:
        raise ValueError(f"Invalid orientation: {orientation}")

# %%
    # We execute seqkit to extract the result.
    # result.stdout   # The command's output (text)
    # result.stderr   # Any error messages
    # result.returncode  # 0 if OK, non-zero if error
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error extracting sequence for edge {edge}: {result.stderr}")
        exit(1)

# %%
    return result.stdout

# %%
# Extract sequences for each circular path and save as concatenated FASTA files.
#
# Use:
# extract_circular_paths_to_fasta(args.gfa, filtered_circular_paths, args.out)
def extract_circular_paths_to_fasta(gfa_file, circular_paths, output_dir):

# %%
    if is_test:
        gfa_file = args.gfa
        circular_paths = filtered_circular_paths
        output_dir = args.out
        print("Test is on.")

# %%
    # Step 1: Convert GFA to FASTA
    fasta_file = os.path.splitext(gfa_file)[0] + ".fa"
    print(f"Converting GFA to FASTA: {fasta_file}")
    run_command(f"gfatools gfa2fa {gfa_file} > {fasta_file}")

# %%
    # Step 2: Process each circular path
    os.makedirs(output_dir, exist_ok=True)
    path_lengths = []
    total_nodes = 0

# %%
    for idx, path in enumerate(circular_paths, start=1):

# %%
        if is_test:
            idx, path = next(enumerate(circular_paths, start=1))
            print(idx)
            print(path)

# %%
        output_fasta = os.path.join(output_dir, f"circular_path_{idx}.fa")
        print(f"Extracting circular path {idx}: {path}")

# %%
        path_sequences = []
        total_nodes += len(path)

# %%
        for edge in path:
            if edge.endswith("+"):
                sequence = extract_edge_sequence(fasta_file, edge, "+")
            elif edge.endswith("-"):
                sequence = extract_edge_sequence(fasta_file, edge, "-")
            else:
                raise ValueError(f"Invalid edge format: {edge}")
            path_sequences.append(sequence.strip())

# %%
        total_length = sum(len(seq) for seq in path_sequences)
        path_lengths.append(total_length)

# %%
        with open(output_fasta, "w") as f:
            for seq in path_sequences:
                f.write(seq + "\n")

# %%
        concatenated_fasta = os.path.join(
            output_dir, f"circular_path_{idx}_concatenated.fa"
        )

# %%
        run_command(f'echo ">path{idx}" >"{concatenated_fasta}"')
        run_command(
            f'seqkit fx2tab {output_fasta} | cut -f2 | tr -d "\n" >>"{concatenated_fasta}"'
        )
        print(f"Circular path {idx} saved as concatenated FASTA: {concatenated_fasta}")

# %%
    average_nodes = round(total_nodes / len(circular_paths), 3) if circular_paths else 0
    circular_path_nodes_file = os.path.join(output_dir, "circular_path_nodes.txt")
    with open(circular_path_nodes_file, "w") as nodes_file:
        nodes_file.write(f"{average_nodes}\n")


# %%
def main():

# %%
    parser = argparse.ArgumentParser(
        description="Extract sequences for circular paths from a GFA file."
    )
    parser.add_argument(
        "-g",
        "--gfa",
        required=True,
        help="Input GFA file")
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Output folder")
    parser.add_argument(
        "-s",
        "--seed",
        required=False,
        help="Seed file containing a list of edges"
    )


# %%
    if is_test:
        # Test arguments
        custom_args = [
            '--gfa', 'input/run-polap-py-find-plastid-gfa2fasta.input1.5.gfa', 
            '--seed', 'input/run-polap-py-find-plastid-gfa2fasta.input2.5.mt.contig.name', 
            '--out', 'output/run-polap-py-find-plastid-gfa2fasta.output.5.52-mtdna'
        ]
        args = parser.parse_args(custom_args)
        print(args)


# %%
    args = parser.parse_args()

# %%
    G = parse_gfa_to_graph(args.gfa)

# %%
    all_circular_paths = find_circular_paths(G)
    print(f"Total number of circular paths found: {len(all_circular_paths)}")

# %%
    filtered_circular_paths = filter_circular_paths_by_seed(
        all_circular_paths, args.seed
    )
    print(f"Number of circular paths after filtering: {len(filtered_circular_paths)}")

# %%
    # Save the count of filtered circular paths to circular_path_count.txt
    circular_path_count_file = os.path.join(args.out, "circular_path_count.txt")
    with open(circular_path_count_file, "w") as count_file:
        count_file.write(f"{len(filtered_circular_paths)}\n")

# %%
    # Exit main if the circular_path_count is too many: i.e., 30
    _fixed_upper_bound_number_segments = 30
    if len(filtered_circular_paths) >= _fixed_upper_bound_number_segments:
        return 1

# %%
    # Save filtered circular paths to circular_path.txt
    circular_path_file = os.path.join(args.out, "circular_path.txt")
    with open(circular_path_file, "w") as path_file:
        for path in filtered_circular_paths:
            path_str = ", ".join(path)
            path_file.write(f"{path_str}\n")

# %%
    # Call the extraction function
    extract_circular_paths_to_fasta(args.gfa, filtered_circular_paths, args.out)

# %%
if __name__ == "__main__":
    main()
