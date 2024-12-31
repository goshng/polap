import os
import subprocess
import sys
import networkx as nx
import argparse


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
    """
    simple_cycles = list(nx.simple_cycles(G))
    circular_components = []

    for cycle in simple_cycles:
        # Include all cycles, including single-node self-loops
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
    for idx, path in enumerate(circular_paths, start=1):
        # Step 2.1: Generate individual FASTA for the path
        output_fasta = os.path.join(output_dir, f"circular_path_{idx}.fa")
        print(f"Extracting circular path {idx}: {path}")

        # Collect sequences for the path
        path_sequences = []
        for edge in path:
            if edge.endswith("+"):
                sequence = extract_edge_sequence(fasta_file, edge, "+")
            elif edge.endswith("-"):
                sequence = extract_edge_sequence(fasta_file, edge, "-")
            else:
                raise ValueError(f"Invalid edge format: {edge}")
            path_sequences.append(sequence.strip())

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
    print(f"Number of circular paths: {len(circular_paths)}")
    for idx, path in enumerate(circular_paths, start=1):
        print(f"Circular Path {idx}: {path}")

    # Call the extraction function
    extract_circular_paths_to_fasta(args.gfa, circular_paths, args.out)


if __name__ == "__main__":
    main()
