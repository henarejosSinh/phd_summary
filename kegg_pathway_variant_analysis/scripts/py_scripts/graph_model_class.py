# ihc.europa@gmail.com
# create a graph using kegg nodes and edges objects
# we'll use networkx 

import networkx as nx
from networkx.readwrite import graph6
import csv

from KeggNode_class import KeggNode
from KeggEdge_class import KeggEdge


class kegg_graph:

    def build_graph(self, file):
        with open(file, "r") as csvfile:
            reader = csv.reader(csvfile)
            next(reader) # skip reader
            for line in reader:
                # get kegg node1-node2 -subtype relationships
                line = line.strip("\n")
                entry1, entry2, _ , subtype = line.split("\t") 

                # check if subtypes in dict, else throw an error
                
                # transform entries to KeggNode obj
                node1 = KeggNode(entry1)
                node2 = KeggNode(entry2)
                G = nx.Graph()
                
                # add nodes to graph
                G.add_node(node1)
                G.add_node(node2)
                
                # model subtype sign
                
                sign = "?"
                if subtype == "activation":
                    sign = "+"
                elif subtype == "inhibition":
                    sign = "-"
                else:
                    sign = "?"
                
                edge = KeggEdge(sign)
                G.add_edge(edge)
                
            return G
        
    @staticmethod
    def parse_subtypes(pathfile):
        subtypes = {}
    
        parser = open(pathfile, "r")
        for line in parser:
            line = line.strip("\n")
            subtypes = {k:v for k,v in line.split("\t")}
            
    @staticmethod
    def build_subpathways(graph):
        startNodes = set()
        endNodes = set()
        allSinglePaths = []
        
        # save start and end nodes of each subpathway
        for i in list(graph.nodes):
            if graph.out_degree(i) == 0 and graph.in_degree(i) >= 1:
                endNodes.add(i)
            if graph.in_degree(i) == 0 and graph.out_degree(i) >= 1:
                startNodes.add(i)
        
        for i in startNodes:
            for j in endNodes:
                allSinglePaths.append(graph.all_single_paths())
        
        return allSinglePaths
    
    # cycles? 
            
        
