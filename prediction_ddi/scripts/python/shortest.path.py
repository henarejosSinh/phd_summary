#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 12:23:33 2020

@author: sinh
"""

import csv
from igraph import *
import argparse

# For functional approach 

'''
Shortest path between two given nodes in a network
'''

#  Parse arguments
parser = argparse.ArgumentParser(description="Process a SIF file to create a "
                                             "graph from which shortest paths "
                                             "will be calculated on the nodes"
                                             " found in the ATT file passed  ")
parser.add_argument("-SIF", metavar="SIF_file", type=str, help="SIF file to "
                                                               "create graph")
parser.add_argument("-ATT", metavar="ATT_file", type=str, help="ATT file with "
                                                               "nodes to "
                                                               "find shortest "
                                                               "paths")
parser.add_argument("-N", metavar="NEIGHBORS", type=int, help="Number of "
                                                              "maximum "
                                                              "neighbors to "
                                                              "add"
                                                              " to paths")
parser.add_argument("-OUT", metavar="OUT_file", type=str, help="File where "
                                                               "modified SIF "
                                                               "will be "
                                                               "written")
args = parser.parse_args()


def get_edges_tuples(edge):  # function to get from original sif directed
    # interaction between
    return g.vs[edge.tuple[0]]['name'], g.vs[edge.tuple[1]]['name']


g = Graph(directed=True)  # works in console
interaction_type = []
pathway = []

with open(args.SIF) as csv_file:
          # newline='') as csv_file:  # network to use for creating an object
          # of graph class

    csv_reader = csv.reader(csv_file, delimiter='\t',
                            quotechar='|')  # must change delimiter depending
          # on case to case
    next(csv_reader)  # skip first line

    for line in csv_reader:
        # nodes
        g.add_vertices(line[0])
        g.add_vertices(line[1])

        # edges
        g.add_edges([(line[0], line[1])])

        # edge attributes
        if len(line) > 1:
            interaction_type.append(line[2])
            pathway.append(line[3])

# edge attributes
g.es['interaction_type'] = interaction_type
g.es['pathway'] = pathway

nodes_to_use = []
dictionary_A_to_B = dict()

with open(args.ATT, newline='') as csv_file:

    csv_reader = csv.reader(csv_file, delimiter=',', quotechar='|')  #
    for line in csv_reader:  # to get index of nodes names column
        # print(line)
        index = line.index('"name"')
        break
    next(csv_reader)
    for line in csv_reader:
        line_str = line[index].replace('"', "")  # Index needs to change
        # appropriately to attribute file passed
        line_str = line_str.replace("'", "")
        nodes_to_use.append(line_str)

for g1 in nodes_to_use:
    for g2 in nodes_to_use:
        if g1 != g2:  # and \

            key = g1 + "_" + g2
            shortest_path = g.get_shortest_paths(g1, to=g2, weights=None,
                                                 mode="out")
            #print(shortest_path)

            if len(shortest_path[0]) <= (args.N + 2) and \
                    len(shortest_path[0]) != 0:
                list_names = []
                for i in shortest_path[0]:  # shortest path uses the
                    # names assigned by Graph
                    name = g.vs[int(i)].attributes()['name']  # return the real
                    # name from the vertex ex 1 = hsa:xx
                    list_names.append(name)  # path from g1 to g2
                dictionary_A_to_B[key] = list_names


# build new sif:

output = open(args.OUT, 'w')
output.write("name1" + "\t" + "name2" + "\t" + "subtype" + "\t" + "pathway")
output.write("\n")

list_node1_to_node2 = []
list_sif = []

for nodes, path in dictionary_A_to_B.items():
    list_node1_to_node2.append(nodes)  # keys
    for j in range(len(path) - 1):
        if [path[j], path[j + 1]] in list_sif or [path[j + 1], path[j]] in \
                list_sif:
            # if A to B or reverse already in sif;
            continue  # continue and don't add to new sif
        list_sif.append([path[j], path[j + 1]])  # else add interactions

        # checkpoint to coerce relations given the data in original sif :
        #print(path[j], path[j + 1])
        source = g.vs.select(name=path[j])  # return given name in graph for
        # the vertex
        target = g.vs.select(name=path[j + 1])
        # Here is were the direction between the 2 nodes being evaluated is
        # corrected to the direction in original sif:
        edges = g.es.select(_between=(source, target))
        for edge in edges:
            # obtain edge attributes
            edge_tuple = get_edges_tuples(edge)
            # print("\n", edge_tuple)
            output.write(
                edge_tuple[0] + "\t" + edge_tuple[1] + "\t" +
                edge["interaction_type"] + "\t" + edge["pathway"])
            output.write("\n")
