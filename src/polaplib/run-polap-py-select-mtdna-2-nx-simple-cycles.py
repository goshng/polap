# Python program to print connected
# components in an undirected graph

# https://stackoverflow.com/questions/54536574/finding-all-circular-path-for-undirected-graph-created-using-networkx
import sys
import pandas as pd
import networkx as nx
import csv

import os
debug = os.getenv("_POLAP_DEBUG", "0")

# Driver Code
if __name__ == "__main__":

    # Create a graph given in the above diagram
    # 5 vertices numbered from 0 to 4

    # g = Graph(5)
    if len(sys.argv) != 3:
        print("Usage: python script.py <path_to_file> <outfile>")
        sys.exit(1)

    file_path = sys.argv[1]
    file_path_out = sys.argv[2]

    # Reading a TSV file into a DataFrame
    df = pd.read_csv(file_path, delimiter="\t")

    # Display the first few rows of the DataFrame
    # print(df.head())

    # df = pd.DataFrame(
    #     {
    #         "Source": [0, 1, 2, 0, 16, 17, 18, 0, 16, 19, 1, 2],
    #         "Target": [1, 2, 3, 3, 17, 18, 19, 16, 19, 3, 17, 18],
    #     }
    # )

    G = nx.from_pandas_edgelist(
        df, source="Source", target="Target", create_using=nx.DiGraph
    )

    # data = nx.cycles.find_cycle(G, orientation="original", source="2+")

    # data = nx.cycles.find_cycle(G, orientation="original")  # , source="4-")
    #
    # # Creating a DataFrame from the data
    # dfo = pd.DataFrame(data, columns=["Source", "Target", "Orientation"])
    #
    # dfo.to_csv(file_path_out, sep="\t", index=False)

    cycles = list(nx.simple_cycles(G))

    # Iterate over each cycle, convert to DataFrame, and save to TSV
    for i, cycle in enumerate(cycles):
        # Prepend 'edge_' to each node in the cycle
        # modified_cycle = [f"edge_{node}" for node in cycle]
        modified_cycle = [f"edge_{node[:-1]}\t{node[-1]}" for node in cycle]

        # Convert the modified cycle to a DataFrame
        df = pd.DataFrame(modified_cycle)

        # Define the full file name with the path and prefix
        file_name = f"{file_path_out}_{i}.tsv"

        # Save the DataFrame to a TSV file without a header and index
        df.to_csv(
            file_name,
            sep="\t",
            index=False,
            header=False,
            quoting=csv.QUOTE_NONE,
            escapechar=" ",
        )

        # print(f"Saved {file_name}")
