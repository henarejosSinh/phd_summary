# Script to detect modules in correlation networks

# from community import generate_dendrogram
# from community import *
import csv
from networkx.algorithms.community import greedy_modularity_communities
import networkx as nx
import pprint as p

import matplotlib.pyplot as plt

G = nx.Graph()

print("Reading File")

# with open("filtered_corr_test_alternative_mse.tsv") as f:
with open("filtered_corr_test_alternative_mse.tsv") as f:
    reader = csv.reader(f, delimiter='\t')
    next(reader, None)  # skip the headers

    for row in reader:

        Gene_id_1 = row[0]
        Gene_id_2 = row[1]

        G.add_edge(Gene_id_1, Gene_id_2)
f.close()

c = list(greedy_modularity_communities(G))

p = nx.shortest_path(G)


C_best = community_louvain.best_partition(G)
dendro = community_louvain.generate_dendrogram(G)
C_0 = community_louvain.partition_at_level(dendro, 0)


M = community_louvain.modularity(C_best, G)
partition = C_best

size = float(len(set(partition.values())))
pos = nx.spring_layout(G)
count = 0.
for com in set(partition.values()) :
    count += 1.
    list_nodes = [nodes for nodes in partition.keys()
                                if partition[nodes] == com]
    nx.draw_networkx_nodes(G, pos, list_nodes, node_size = 20,
                                node_color = str(count / size))
nx.draw_networkx_edges(G, pos, alpha=0.5)
plt.show()

w = open("clusters_mse.tsv", "w")
n = 0
for i in c:
    n += 1
    for element in i:
        print(element)
        w.write("\t{}\t{}\n".format(str(n), element))