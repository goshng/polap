# %%
msg = "Hello World"
print(msg)

# %%
msg = "Hello again"
print(msg)

# %%
import networkx as nx
import os
G = nx.Graph()
all_nodes = set()

# %%
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
# %%
gfa = "1.gfa"
graph, all_nodes = parse_gfa_to_graph(gfa)
# %%
print(graph)
print(all_nodes)
# %%
file_object = open("2.txt", "r")
content = file_object.readline()
print(content)

# %%
content = file_object.readline()
print(content)
# %%
file_object.close()
#%%
with open("2.txt", "r") as file_object:
    content = file_object.read()
    print(content)
# File is automatically closed here
#%%
