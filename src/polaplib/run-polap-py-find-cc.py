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

# Python script to print connected components in an undirected graph
#
# Links:
# 1. IO processing
# https://chatgpt.com/c/66dbce52-c6ec-800e-a940-629c6a333d31
#
# 2. DFS(Depth First Search)
# https://www.geeksforgeeks.org/python-program-for-depth-first-search-or-dfs-for-a-graph/

import sys


class Graph:

    # init function to declare class variables
    def __init__(self, V):
        self.V = V
        self.adj = [[] for i in range(V)]

    def DFSUtil(self, temp, v, visited):

        # Mark the current vertex as visited
        visited[v] = True

        # Store the vertex to list
        temp.append(v)

        # Repeat for all vertices adjacent
        # to this vertex v
        for i in self.adj[v]:
            if visited[i] == False:

                # Update the list
                temp = self.DFSUtil(temp, i, visited)
        return temp

    # method to add an undirected edge
    def addEdge(self, v, w):
        self.adj[v].append(w)
        self.adj[w].append(v)

    # Method to retrieve connected components
    # in an undirected graph
    def connectedComponents(self):
        visited = []
        cc = []
        for i in range(self.V):
            visited.append(False)
        for v in range(self.V):
            if visited[v] == False:
                temp = []
                cc.append(self.DFSUtil(temp, v, visited))
        return cc


def read_search_numbers(file_path):
    # Read numbers from the text file
    search_numbers = []
    with open(file_path, "r") as file:
        for line in file:
            # Convert the line to an integer and add to the list
            search_numbers.append(int(line.strip()))
    return search_numbers


def find_arrays_with_numbers(array_of_arrays, search_numbers):
    # Find all sub-arrays that contain any of the search numbers
    result_arrays = []
    for sub_array in array_of_arrays:
        if any(num in sub_array for num in search_numbers):
            result_arrays.append(sub_array)
    return result_arrays


def flatten_and_unique(arrays):
    # Flatten the list of lists and remove duplicates by converting to a set
    flattened = [item for sublist in arrays for item in sublist]
    unique_elements = set(flattened)
    return unique_elements


def write_unique_elements_to_file(unique_elements, file_path):
    # Write unique elements to a text file, one per line
    with open(file_path, "w") as file:
        for element in sorted(unique_elements):  # Sort for consistent output
            file.write(f"{element}\n")


# def addEdge(node1, node2):
#     # Example function definition (you can replace this with your actual implementation)
#     print(f"Edge added between {node1} and {node2}")


# Replace 'edges.txt' with the path to your text file
# process_edges("edges.txt")

# Driver Code
if __name__ == "__main__":

    # Create a graph given in the above diagram
    # 5 vertices numbered from 0 to 4

    # g = Graph(5)
    if len(sys.argv) != 4:
        print("Usage: python script.py <path_to_file> <search_numbers_file> <outfile>")
        sys.exit(1)

    file_path = sys.argv[1]
    file_path_search_number = sys.argv[2]
    file_path_out = sys.argv[3]

    unique_numbers = set()

    with open(file_path, "r") as file:
        for line in file:
            numbers = line.strip().split()
            if len(numbers) == 2:
                node1, node2 = map(int, numbers)
                unique_numbers.update([node1, node2])  # Add the numbers to the set

    # Call the Graph function with the number of unique nodes
    g = Graph(len(unique_numbers))

    with open(file_path, "r") as file:
        for line in file:
            # Split the line into two numbers, assuming they are space-separated
            numbers = line.strip().split()
            if len(numbers) == 2:
                node1, node2 = map(int, numbers)  # Convert the two numbers to integers
                g.addEdge(node1, node2)
            else:
                print(f"Skipping invalid line: {line}")

    # g.addEdge(1, 0)
    # g.addEdge(2, 1)
    # g.addEdge(3, 4)
    # not using networkx
    cc = g.connectedComponents()
    # print("Following are connected components")
    # print(cc)

    search_numbers = read_search_numbers(file_path_search_number)

    # Find arrays containing any of the search numbers
    result_arrays = find_arrays_with_numbers(cc, search_numbers)

    # Flatten the resulting arrays and get unique elements
    unique_elements = flatten_and_unique(result_arrays)

    # print("Unique elements from the arrays containing any of the search numbers:")
    # print(unique_elements)

    write_unique_elements_to_file(unique_elements, file_path_out)
