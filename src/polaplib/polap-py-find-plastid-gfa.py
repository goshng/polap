#!/usr/bin/env python3

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

import argparse
import os
import sys
import networkx as nx

# Assumes these helper functions are defined:
# - parse_gfa_with_inferred_links
# - find_connected_component
# - is_plastid_path
# - extract_plastid_sequence
# - print_component_edges_by_role
# - write_fasta


def write_edge_list(path_with_orientations, outfile):
    unique_edges = sorted({name for name, _ in path_with_orientations})
    with open(outfile, "w") as f:
        for name in unique_edges:
            f.write(name + "\n")


def write_fasta_to_stream(sequence, output_path, seq_id="plastid"):
    if output_path == "-":
        out = sys.stdout
    else:
        out = open(output_path, "w")

    with out:
        out.write(f">{seq_id}\n")
        for i in range(0, len(sequence), 70):
            out.write(sequence[i : i + 70] + "\n")


def reverse_orient(orient):
    return "+" if orient == "-" else "-"


def reverse_link(from_seg, from_orient, to_seg, to_orient):
    return (
        f"{to_seg}{reverse_orient(to_orient)}",
        f"{from_seg}{reverse_orient(from_orient)}",
    )


def parse_gfa_with_inferred_links(gfa_path):
    segments = {}
    links = []
    with open(gfa_path) as f:
        for line in f:
            if line.startswith("S"):
                parts = line.strip().split("\t")
                segments[parts[1]] = parts[2]
            elif line.startswith("L"):
                parts = line.strip().split("\t")
                from_seg, from_orient = parts[1], parts[2]
                to_seg, to_orient = parts[3], parts[4]
                # Explicit link
                a = f"{from_seg}{from_orient}"
                b = f"{to_seg}{to_orient}"
                links.append((a, b))
                # Inferred reverse link
                rev_a, rev_b = reverse_link(from_seg, from_orient, to_seg, to_orient)
                links.append((rev_a, rev_b))
    return segments, links


def find_connected_component(gfa_path, start_edge):
    segments, links = parse_gfa_with_inferred_links(gfa_path)
    G = nx.Graph()
    G.add_edges_from(links)

    if start_edge not in G:
        print(f"Start edge {start_edge} not found.")
        return segments, None, []

    for component in nx.connected_components(G):
        if start_edge in component:
            print("Connected component edges:")
            for node in sorted(component):
                print(f"  {node}")
            return segments, G.subgraph(component).copy(), sorted(component)

    return segments, None, []


def is_plastid_path(component_graph, segments, verbose=False):
    """
    Detect plastid genome paths:
    - Type 1: Single-edge circular: self-loop on + only (edgeX+ -> edgeX+)
    - Type 2: Tripartite: LSC -> IR+ -> SSC -> IR- (3 distinct edges, IR bidirectional, with verified connectivity)
    """
    from collections import defaultdict

    nodes = set(component_graph.nodes)

    # Build orientation map
    edge_orients = defaultdict(set)
    for node in nodes:
        if node[-1] in "+-":
            edge_orients[node[:-1]].add(node[-1])

    if verbose:
        print("\n[Edge orientation map]")
        for name, orients in edge_orients.items():
            print(f"  {name}: {sorted(orients)}")

    # ---- Type 1: Single-edge circular (self-loop on + only) ----
    for name in edge_orients:
        plus = f"{name}+"
        if plus in component_graph and component_graph.has_edge(plus, plus):
            if verbose:
                print(
                    f"\n[Detected single-edge circular plastid using {plus} self-loop]"
                )
            return [(name, "+")]

    # ---- Type 2: LSC–IR+–SSC–IR− ----
    candidate_irs = [name for name, o in edge_orients.items() if "+" in o and "-" in o]

    if verbose:
        print(f"\n[Candidate IRs]: {candidate_irs}")

    for ir in candidate_irs:
        ir_plus = f"{ir}+"
        ir_minus = f"{ir}-"

        neighbors_plus = set(component_graph.neighbors(ir_plus))
        neighbors_minus = set(component_graph.neighbors(ir_minus))

        neighbors_plus.discard(ir_minus)
        neighbors_minus.discard(ir_plus)

        if verbose:
            print(f"\nIR+ neighbors of {ir_plus}: {neighbors_plus}")
            print(f"IR- neighbors of {ir_minus}: {neighbors_minus}")

        for a_node in neighbors_plus:
            for b_node in neighbors_minus:
                a_name, a_orient = a_node[:-1], a_node[-1]
                b_name, b_orient = b_node[:-1], b_node[-1]

                if a_name in (ir,) or b_name in (ir,) or a_name == b_name:
                    continue
                if a_name not in segments or b_name not in segments:
                    continue

                len_a, len_b = len(segments[a_name]), len(segments[b_name])
                lsc, ssc = (
                    ((a_name, a_orient), (b_name, b_orient))
                    if len_a >= len_b
                    else ((b_name, b_orient), (a_name, a_orient))
                )

                if verbose:
                    print("\n[Detected 3-edge plastid path with roles]:")
                    print(f"  LSC:  {lsc[0]}{lsc[1]}")
                    print(f"  IR+:  {ir}+")
                    print(f"  SSC:  {ssc[0]}{ssc[1]}")
                    print(f"  IR-:  {ir}-")

                return [lsc, (ir, "+"), ssc, (ir, "-")]

    if verbose:
        print("\n[No valid plastid structure found]")
    return None


def reverse_complement(seq):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]


def extract_plastid_sequence(segments, path_with_orientations):
    sequence = ""
    for seg, orient in path_with_orientations:
        seq = segments[seg]
        sequence += seq if orient == "+" else reverse_complement(seq)
    return sequence


def print_component_edges_by_role(path_with_orientations):
    """
    Print unique edge names in order:
    - Single-edge case: label as Circular
    - Tripartite case: LSC, IR, SSC
    """
    if not path_with_orientations:
        print("No valid path found.")
        return

    unique_names = []
    seen = set()
    for name, _ in path_with_orientations:
        if name not in seen:
            unique_names.append(name)
            seen.add(name)

    print("Edges in plastid component (by role):")
    if len(unique_names) == 1:
        print(f"  Circular plastid edge: {unique_names[0]}")
    elif len(unique_names) == 3:
        print(f"  LSC:  {unique_names[0]}")
        print(f"  IR:   {unique_names[1]}")
        print(f"  SSC:  {unique_names[2]}")
    else:
        for i, name in enumerate(unique_names):
            print(f"  Region{i+1}: {name}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract plastid genome path and sequence from GFA."
    )
    parser.add_argument("gfa", help="Input GFA file")
    parser.add_argument(
        "--seed", default="edge_1", help="Seed edge name (default: edge_1)"
    )
    parser.add_argument(
        "--mtcontig",
        default="mt.contig.name-1",
        help="Output file for plastid edge list",
    )
    parser.add_argument(
        "--fasta",
        default="-",
        help="Output FASTA file for plastid sequence (default: stdout)",
    )

    args = parser.parse_args()

    if not os.path.isfile(args.gfa):
        parser.error(f"GFA file not found: {args.gfa}")

    # Add '+' if no orientation provided in --seed
    seed = args.seed if args.seed[-1] in "+-" else args.seed + "+"

    # Step 1: Parse and find component
    segments, graph, component = find_connected_component(args.gfa, seed)
    if not graph:
        print("No connected component found.")
        return

    # Step 2: Detect plastid path
    path = is_plastid_path(graph, segments, verbose=True)
    if not path:
        print("No valid plastid path found.")
        return

    # Step 3: Output results
    print_component_edges_by_role(path)

    if args.mtcontig:
        write_edge_list(path, args.mtcontig)
        print(f"Wrote plastid edge list to {args.mtcontig}")

    if args.fasta:
        seq = extract_plastid_sequence(segments, path)
        write_fasta_to_stream(seq, args.fasta)
        if args.fasta != "-":
            print(f"Wrote plastid FASTA to {args.fasta}")


if __name__ == "__main__":
    main()
