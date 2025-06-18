# polaplib/run-polap-py-find-cc-with-seeds.py
#
# TODO: document
#
# This script parses gfa to extract plastid DNA.
#
# Used by:
# function _run_polap_step-disassemble-seeds-graph {

import argparse
import networkx as nx

import os
debug = os.getenv("_POLAP_DEBUG", "0")

def parse_gfa_to_graph(gfa_file):
    """
    Parse the GFA file and construct an undirected graph.
    Also collect all nodes from lines starting with "S".
    """
    G = nx.Graph()
    all_nodes = set()

    with open(gfa_file, "r") as file:
        for line in file:
            if line.startswith("S"):
                fields = line.strip().split("\t")
                node = fields[1]
                all_nodes.add(node)
            elif line.startswith("L"):
                # Parse link lines
                fields = line.strip().split("\t")
                if len(fields) >= 3:
                    source = fields[1]
                    target = fields[3]
                    G.add_edge(source, target)

    return G, all_nodes


def load_nodes_from_file(node_file):
    """
    Load the list of nodes from the given file. Extract only the first item (edge_<number>) from each line.
    """
    with open(node_file, "r") as file:
        nodes = {line.strip().split("\t")[0] for line in file if line.strip()}
    return nodes


def find_connected_components(graph, nodes_of_interest):
    """
    Find all connected components that contain at least one of the nodes of interest.
    Include nodes from the nodes_of_interest as single-node components if they are isolated.
    """
    connected_components = []
    processed_nodes = set()

    # Find components containing nodes of interest
    for component in nx.connected_components(graph):
        if (
            nodes_of_interest & component
        ):  # Check for intersection with nodes_of_interest
            connected_components.append(component)
            processed_nodes.update(component)

    # Include nodes not part of any component
    isolated_nodes = nodes_of_interest - processed_nodes
    for node in isolated_nodes:
        connected_components.append({node})  # Add as a single-node component

    return connected_components


def main():
    parser = argparse.ArgumentParser(
        description="Find connected components containing specific nodes in a GFA graph."
    )
    parser.add_argument("--gfa", required=True, help="Input GFA file.")
    parser.add_argument(
        "--nodes",
        required=True,
        help="File containing the list of nodes (one per line).",
    )
    parser.add_argument(
        "--output",
        required=False,
        help="Output file to save the results.",
        default="connected_components.txt",
    )

    args = parser.parse_args()

    # Parse the GFA file and construct the graph
    print("Parsing GFA file...")
    graph, all_nodes = parse_gfa_to_graph(args.gfa)

    # Load the list of nodes
    print("Loading nodes of interest...")
    nodes_of_interest = load_nodes_from_file(args.nodes)

    # Validate that all nodes of interest are part of the graph
    missing_nodes = nodes_of_interest - all_nodes
    if missing_nodes:
        print(
            f"Error: The following nodes are not in the graph: {', '.join(missing_nodes)}"
        )
        exit(1)

    # Find connected components containing at least one node of interest
    print("Finding connected components...")
    connected_components = find_connected_components(graph, nodes_of_interest)

    # Save or display the results
    if connected_components:
        print(f"Found {len(connected_components)} connected components.")
        with open(args.output, "w") as output_file:
            for idx, component in enumerate(connected_components, start=1):
                component_str = ", ".join(component)
                output_file.write(f"{component_str}\n")
        print(f"Results saved to {args.output}")
    else:
        print("No connected components found containing the specified nodes.")


if __name__ == "__main__":
    main()
