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

# polaplib/polap-py-find-cc-with-seeds.py

################################################################################
# This script parses gfa to extract plastid DNA edges.
# It uses networkx python library to find connected components.
# We have gfa and a list of edges that the connected components will have. 
# Because networkx is a general graph library, we need to convert the gfa
# to a more general graph form or edges and nodes.
#
# Files:
# ======
#
# input/run-polap-py-find-cc-with-seeds.input1.2.gfa
# --------------------------------------------------
# S       edge_1  *       LN:i:155591     dp:i:27
# L       edge_1  +       edge_1  +       0M      L1:i:155591     L2:i:155591     RC:i:0
#
# input/run-polap-py-find-cc-with-seeds.input2.2.txt
# --------------------------------------------------
# edge_1  155591  27      1       32      56      1       2778
#
# output/run-polap-py-find-cc-with-seeds.output1.2.txt
# ----------------------------------------------------
# edge_1
#
# Example:
# python polap-py-find-cc-with-seeds.py \
#   --gfa input/run-polap-py-find-cc-with-seeds.input1.2.gfa \
#   --nodes input/run-polap-py-find-cc-with-seeds.input2.2.txt  \
#   --output output/run-polap-py-find-cc-with-seeds.output1.2.txt
#
# Additional Files:
# =================
#
# input/run-polap-py-find-cc-with-seeds.input1.5.gfa
# --------------------------------------------------
# S       edge_1  *       LN:i:86832      dp:i:58
# S       edge_2  *       LN:i:17970      dp:i:66
# S       edge_4  *       LN:i:25975      dp:i:128
# L       edge_1  +       edge_4  +       0M      L1:i:86832      L2:i:25975      RC:i:59
# L       edge_1  -       edge_4  +       0M      L1:i:86832      L2:i:25975      RC:i:55
# L       edge_2  +       edge_4  -       0M      L1:i:17970      L2:i:25975      RC:i:69
# L       edge_2  -       edge_4  -       0M      L1:i:17970      L2:i:25975      RC:i:46
#
# output/run-polap-py-find-cc-with-seeds.output1.5.txt
# ----------------------------------------------------
# edge_2, edge_4, edge_1
#
# Used by:
# function _run_polap_step-disassemble-seeds-graph {
#
################################################################################

# %%
import argparse
import networkx as nx
import os

debug = os.getenv("_POLAP_DEBUG", "0")

# %%
is_test = True
is_test = False

# %%

# Convert gfa to a graph form so that we can use networkx module to find
# the connected component.
# Parse the GFA file and construct an undirected graph.
# Also collect all nodes from lines starting with "S".
#
# Use:
# graph, all_nodes = parse_gfa_to_graph(args.gfa)
def parse_gfa_to_graph(gfa_file):

## %%
    if is_test:
        gfa_file = args.gfa
        print("Test is on.")

## %%
    # An empty graph: G
    G = nx.Graph()
    # An empty python set: all_nodes
    all_nodes = set()

## %%
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

## %%
    return G, all_nodes

# %%
# Load the list of nodes from the given file. 
# Extract only the first item (edge_<number>) from each line.
def load_nodes_from_file(node_file):

## %%
    if is_test:
        node_file = args.nodes
        print("Test is on.")
## %%
    with open(node_file, "r") as file:
        nodes = {line.strip().split("\t")[0] for line in file if line.strip()}
## %%
    return nodes

# %%

# Find all connected components that contain at least one of the nodes of interest.
# Include nodes from the nodes_of_interest as single-node components if they are isolated.
#
# Use:
# connected_components = find_connected_components(graph, nodes_of_interest)
def find_connected_components(graph, nodes_of_interest):

## %%
    if is_test:
        print("Test is on.")

## %%
    connected_components = []
    processed_nodes = set()

## %%
    # Find components containing nodes of interest
    for component in nx.connected_components(graph):
        if (
            nodes_of_interest & component
        ):  # Check for intersection with nodes_of_interest
            connected_components.append(component)
            processed_nodes.update(component)

## %%
    # Include nodes not part of any component
    isolated_nodes = nodes_of_interest - processed_nodes
    for node in isolated_nodes:
        connected_components.append({node})  # Add as a single-node component

## %%
    return connected_components


# %%
def main():

# %%
    parser = argparse.ArgumentParser(
        description="Find connected components containing specific nodes in a GFA graph."
    )
    parser.add_argument(
        "--gfa",
        required=True,
        help="Input GFA file."
    )
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

# %%
    if is_test:
        # Test arguments
        custom_args = [
            '--gfa', 'input/run-polap-py-find-cc-with-seeds.input1.2.gfa',
            '--nodes', 'input/run-polap-py-find-cc-with-seeds.input2.2.txt', 
            '--output', 'output/run-polap-py-find-cc-with-seeds.output1.2.txt'
        ]
        args = parser.parse_args(custom_args)
        print(args)

# %%
    args = parser.parse_args()

# %%
    # Parse the GFA file
    graph, all_nodes = parse_gfa_to_graph(args.gfa)

# %%
    # Load the list of nodes for finding connected components
    nodes_of_interest = load_nodes_from_file(args.nodes)

# %%
    # Validate that all nodes of interest are part of the graph
    missing_nodes = nodes_of_interest - all_nodes
    if missing_nodes:
        print(
            f"Error: The following nodes are not in the graph: {', '.join(missing_nodes)}"
        )
        exit(1)

# %%
    # Find connected components containing at least one node of interest
    connected_components = find_connected_components(graph, nodes_of_interest)

# %%
    # Save or display the results
    if connected_components:
        with open(args.output, "w") as output_file:
            for idx, component in enumerate(connected_components, start=1):
                component_str = ", ".join(component)
                output_file.write(f"{component_str}\n")
        print(f"Connected components: {args.output}")
    else:
        print("No connected components found containing the specified nodes.")

# %%

if __name__ == "__main__":
    main()
