# polaplib/run-polap-py-find-plastid-gfa2fasta-filtering-edges.py
#
# Used by:
# nothing

import os
import subprocess
import sys
import networkx as nx
import argparse

debug = os.getenv("_POLAP_DEBUG", "0")

def parse_gfa_to_graph_with_reverses_and_filter(gfa_file, seed_file):
    """
    Parse a GFA file and construct a directed graph, generating reverse complement edges.
    Filter the link lines based on a seed file containing edge_<number>.
    """
    G = nx.DiGraph()

    # Load the seed edges
    seed_edges = set()
    if seed_file:
        with open(seed_file, "r") as seed:
            seed_edges = {line.strip() for line in seed if line.strip()}

    with open(gfa_file, "r") as file:
        for line in file:
            if line.startswith("L"):  # Process Link lines
                fields = line.strip().split("\t")
                if len(fields) >= 5:
                    source, source_orient, target, target_orient = fields[1:5]

                    # Filter by seed edges
                    if source in seed_edges or target in seed_edges:
                        # Add the original edge
                        G.add_edge(
                            f"{source}{source_orient}", f"{target}{target_orient}"
                        )

                        # Generate and add the reverse complement edge
                        rev_source_orient = "-" if source_orient == "+" else "+"
                        rev_target_orient = "-" if target_orient == "+" else "+"
                        G.add_edge(
                            f"{target}{rev_target_orient}",
                            f"{source}{rev_source_orient}",
                        )
    return G


def find_circular_paths(G):
    """
    Find all circular paths, including single-node cycles (self-loops).
    Exclude paths with more than 4 nodes.
    """
    simple_cycles = list(nx.simple_cycles(G))
    circular_components = []

    for cycle in simple_cycles:
        if len(cycle) <= 4:  # Exclude paths with more than 4 nodes
            # Extract edge names without orientation (+/-)
            unique_edges = {edge.rstrip("+-") for edge in cycle}
            if len(unique_edges) in [1, 3]:  # Keep only paths with 1 or 3 unique edges
                circular_components.append(cycle)

    return circular_components


def run_command(command):
    """Run a shell command and handle errors."""
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
        exit(1)


def extract_edge_sequence(fasta_file, edge, orientation):
    """
    Extract a sequence for a given edge with the correct orientation.
    """
    base_edge = edge.split("+")[0].split("-")[0]  # Remove signs to get base edge name
    if orientation == "+":
        # Extract sequence without reverse complement
        command = f"seqkit grep -p {base_edge} {fasta_file} | seqkit seq -t dna -v"
    elif orientation == "-":
        # Extract sequence with reverse complement
        command = (
            f"seqkit grep -p {base_edge} {fasta_file} | seqkit seq -t dna -v -r -p"
        )
    else:
        raise ValueError(f"Invalid orientation: {orientation}")

    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error extracting sequence for edge {edge}: {result.stderr}")
        exit(1)
    return result.stdout


def extract_circular_paths_to_fasta(gfa_file, circular_paths, output_dir):
    """
    Extract sequences for each circular path and save as concatenated FASTA files.
    """
    # Step 1: Convert GFA to FASTA
    fasta_file = os.path.splitext(gfa_file)[0] + ".fa"
    print(f"Converting GFA to FASTA: {fasta_file}")
    run_command(f"gfatools gfa2fa {gfa_file} > {fasta_file}")

    # Step 2: Process each circular path
    os.makedirs(output_dir, exist_ok=True)
    path_lengths = []
    total_nodes = 0

    for idx, path in enumerate(circular_paths, start=1):
        # Step 2.1: Generate individual FASTA for the path
        output_fasta = os.path.join(output_dir, f"circular_path_{idx}.fa")
        print(f"Extracting circular path {idx}: {path}")

        # Collect sequences for the path
        path_sequences = []
        total_nodes += len(path)
        for edge in path:
            if edge.endswith("+"):
                sequence = extract_edge_sequence(fasta_file, edge, "+")
            elif edge.endswith("-"):
                sequence = extract_edge_sequence(fasta_file, edge, "-")
            else:
                raise ValueError(f"Invalid edge format: {edge}")
            path_sequences.append(sequence.strip())

        # Calculate the total length of the concatenated sequence
        total_length = sum(len(seq) for seq in path_sequences)
        path_lengths.append(total_length)

        # Write the individual FASTA file
        with open(output_fasta, "w") as f:
            for seq in path_sequences:
                f.write(seq + "\n")

        # Step 2.2: Concatenate the sequences into a single FASTA file
        concatenated_fasta = os.path.join(
            output_dir, f"circular_path_{idx}_concatenated.fa"
        )
        run_command(f'echo ">path{idx}" >"{concatenated_fasta}"')
        run_command(
            f'seqkit fx2tab {output_fasta} | cut -f2 | tr -d "\n" >>"{concatenated_fasta}"'
        )
        print(f"Circular path {idx} saved as concatenated FASTA: {concatenated_fasta}")

    # Save average number of nodes per circular path to a file
    average_nodes = total_nodes / len(circular_paths) if circular_paths else 0
    circular_path_nodes_file = os.path.join(output_dir, "circular_path_nodes.txt")
    with open(circular_path_nodes_file, "w") as nodes_file:
        nodes_file.write(f"{average_nodes}\n")


def main():

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Extract sequences for circular paths from a GFA file."
    )
    parser.add_argument("-g", "--gfa", required=True, help="Input GFA file")
    parser.add_argument("-o", "--out", required=True, help="Output folder")
    parser.add_argument(
        "-s", "--seed", required=False, help="Seed file containing a list of edges"
    )

    args = parser.parse_args()

    # Ensure required arguments are provided
    if not args.gfa or not args.out:
        parser.error("Both --gfa and --out options are required.")

    # Parse the graph with optional filtering by seed file
    G = parse_gfa_to_graph_with_reverses_and_filter(args.gfa, args.seed)

    # List all circular paths (including self-loops)
    circular_paths = find_circular_paths(G)

    # Save circular path count and paths to text files
    circular_path_count_file = os.path.join(args.out, "circular_path_count.txt")
    circular_path_file = os.path.join(args.out, "circular_path.txt")

    with open(circular_path_count_file, "w") as count_file:
        count_file.write(f"{len(circular_paths)}\n")

    with open(circular_path_file, "w") as path_file:
        for idx, path in enumerate(circular_paths, start=1):
            path_str = ", ".join(path)
            path_file.write(f"Circular Path {idx}: {path_str}\n")

    print(f"Number of circular paths: {len(circular_paths)}")
    print(f"Circular paths saved to {circular_path_file}")

    # Call the extraction function
    extract_circular_paths_to_fasta(args.gfa, circular_paths, args.out)


if __name__ == "__main__":
    main()
