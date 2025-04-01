# Python program to find a mtDNA sequence from a GFA file.
#
# https://chatgpt.com/share/672f11cd-782c-800e-b220-615f6b593aa6
#
# ### 1. **Allowing Node Revisit in Cycle Calculation**
# If revisiting nodes is allowed, you can modify your search approach by not requiring that each node is only visited once. This can be handled using the following methods:
#
# #### Modified BFS or Dijkstra’s for Shortest Path Cycles
# 1. **Start from Each Node**:
#    - Use Dijkstra’s algorithm (for weighted graphs) or BFS (for unweighted graphs) to explore paths starting from each node.
#    - Track paths and revisit nodes if they lead to shorter cycles or cycles covering more nodes.
# 2. **Cycle Formation**:
#    - For each path that returns to the starting node, calculate its length and the number of unique nodes it visits.
#    - Keep track of the cycle that has the shortest path length while covering the maximum number of unique nodes.
# 3. **Implementation**:
#    - This method uses Dijkstra’s or BFS but allows revisiting nodes as long as it leads to shorter or more comprehensive cycles.
#
# ### 2. **Using NetworkX to Find Cycles with Node Revisits**
# You can adapt `networkx.simple_cycles` or other path-finding methods to construct cycles that revisit nodes if necessary.
#
# Here’s a sample code for this approach:
#
# ```python
# import networkx as nx
# from networkx.algorithms.shortest_paths.weighted import single_source_dijkstra
#
# def shortest_cycle_max_nodes_revisitable(graph):
#     max_nodes_visited = 0
#     shortest_cycle_length = float('inf')
#     best_cycle = None
#
#     for start_node in graph.nodes:
#         # Using Dijkstra’s or BFS to explore paths from each node
#         visited_paths = {}  # Track shortest paths to revisit nodes
#         queue = [(0, start_node, [start_node])]  # (path length, current node, path)
#
#         while queue:
#             path_length, node, path = queue.pop(0)
#
#             # Form a cycle if possible and check for unique node count
#             if node == start_node and len(path) > 1:
#                 unique_nodes = len(set(path))
#                 if unique_nodes > max_nodes_visited or (unique_nodes == max_nodes_visited and path_length < shortest_cycle_length):
#                     max_nodes_visited = unique_nodes
#                     shortest_cycle_length = path_length
#                     best_cycle = path[:]
#
#             # Explore neighbors, allowing revisits
#             for neighbor in graph.neighbors(node):
#                 new_path_length = path_length + graph[node][neighbor].get('weight', 1)
#                 new_path = path + [neighbor]
#
#                 # Only enqueue paths that improve visited paths or lead to shorter cycles
#                 if neighbor not in visited_paths or new_path_length < visited_paths[neighbor]:
#                     visited_paths[neighbor] = new_path_length
#                     queue.append((new_path_length, neighbor, new_path))
#
#     return best_cycle, shortest_cycle_length
# ```
#
# ### Explanation
# - **Cycle Detection**: Each time the algorithm returns to the `start_node` and forms a cycle, it checks if this cycle is the shortest one with the most unique nodes.
# - **Revisits Allowed**: Paths can revisit nodes to explore new potential cycles that might be shorter or visit more nodes.
# - **Tracking Paths**: A `visited_paths` dictionary ensures only the shortest paths to each node are explored further to avoid redundant computations.
#
# ### Notes
# - **Performance**: Allowing revisits generally increases the search space, so this method is computationally intensive for large graphs.
# - **Cycle Selection**: The algorithm always updates the best cycle if it finds one that is either shorter or has more unique nodes.
#
# This approach leverages node revisits to expand the search, potentially finding cycles that cover more nodes or result in shorter paths.

# https://stackoverflow.com/questions/54536574/finding-all-circular-path-for-undirected-graph-created-using-networkx
import sys
import pandas as pd
import networkx as nx
import csv
from networkx.algorithms.shortest_paths.weighted import single_source_dijkstra


def shortest_cycle_max_nodes_revisitable(graph):
    max_nodes_visited = 0
    shortest_cycle_length = float("inf")
    best_cycle = None

    for start_node in graph.nodes:
        # Using Dijkstra’s or BFS to explore paths from each node
        visited_paths = {}  # Track shortest paths to revisit nodes
        queue = [(0, start_node, [start_node])]  # (path length, current node, path)

        while queue:
            path_length, node, path = queue.pop(0)

            # Form a cycle if possible and check for unique node count
            if node == start_node and len(path) > 1:
                unique_nodes = len(set(path))
                if unique_nodes > max_nodes_visited or (
                    unique_nodes == max_nodes_visited
                    and path_length < shortest_cycle_length
                ):
                    max_nodes_visited = unique_nodes
                    shortest_cycle_length = path_length
                    best_cycle = path[:]

            # Explore neighbors, allowing revisits
            for neighbor in graph.neighbors(node):
                new_path_length = path_length + graph[node][neighbor].get("weight", 1)
                new_path = path + [neighbor]

                # Only enqueue paths that improve visited paths or lead to shorter cycles
                if (
                    neighbor not in visited_paths
                    or new_path_length < visited_paths[neighbor]
                ):
                    visited_paths[neighbor] = new_path_length
                    queue.append((new_path_length, neighbor, new_path))

    return best_cycle, shortest_cycle_length


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
