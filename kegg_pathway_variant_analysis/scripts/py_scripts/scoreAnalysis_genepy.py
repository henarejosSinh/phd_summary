# ihc.europa@gmail.com
# script to calculate the score of a given node

# imports
from pprint import pprint
import statistics
import re
from functools import reduce
from abc import ABC, abstractclassmethod, abstractmethod
import typing

# classes from the project
import EffectScore
import AttEntry
import Path
import Subpathway
from functions import *

# vcf, variant classes
# sys.path.insert(1, '/home/ihenarejos/workspace/projects/pof/scripts/python')
from class_vcf import *
from class_variant import *
from scoreAnalysis import scoreAnalysis


class genepy_score(scoreAnalysis):
    
    def __init__(self, hsa_genes_dict, variants_x_genes, effect_scores, genepy_files):
        super().__init__(hsa_genes_dict, variants_x_genes, effect_scores)
        
        self.genepy_scores = genepy_files
        
    def __get_gene_value(self, gene: str, sample: str) -> typing.Union[float, str]:
        """[summary]

        Args:
            gene (str): [gene to search for a genepy score]
            sample (str): [patient to search for]

        Returns:
            float | str: [gene score]
        """        
        matching = [s for s in self.genepy_scores if sample in s]
        if matching:
            sample_file = open(matching[0], "r")
            gene_dict = {k:v for k,v in (line.strip().split("\t") for line in sample_file)}
            
            gene_val = 1.0            
            if gene in gene_dict.keys():
                gene_val = gene_dict[gene]
    
            print(f'\ngene "{gene}" score: {gene_val}') 
            return gene_val
        else:
            print(f'\nERROR {sample} not found in list of files')
            traceback.print_exc()
    
    def __calculate_node_value(self, node, att_dict,
                                sample):
        formula = None
        if node not in att_dict.keys():
            raise ValueError(f'Node to calculate value: {node} not found in Att entries dictionary')
        
        # use AttEntry class to retrieve the gene names and them find
        # the pattern to extract them and save them in a list
        attributes = att_dict[node]
        formula = attributes.get_gene_names()
        # matcher = self.gene_pattern.search(formula)
        matcher = self.gene_pattern.findall(formula)
        
        genes = []
        # while matcher is not None:
        for gene in matcher:
            genes.append(gene[0])
        
        print(f'''
             .att entry: {attributes.to_string()}
              genes in entry: {formula}
              matcher: {matcher}
              genes: {genes}
              ''')
        
        for gene in genes:
            gene_value = self.__get_gene_value(gene, sample)
            formula = formula.replace(gene, str(gene_value), 1)
        
        return self.eval(formula)  
    
    def apply(self, elem):
        elem : Subpathway
        
        hsa = elem.get_hsa()
        subpathway = elem.get_subpathway()
        sample = elem.get_sample()
        
        # create dict, list to save results
        node_values_dict = {}
        path_results = []
        
        
        for path in subpathway.get_paths():
            score = 1.0
            path_nodes = path.get_path()
            node_values = []
            
            for i in range(0, len(path_nodes)):
                node = path_nodes[i]
                
                if node not in node_values_dict.keys():
                    node_value = self.__calculate_node_value(node, 
                                                            subpathway.get_att_entries(), 
                                                            sample)
                    node_values_dict[node] = node_value
                else:
                    node_value = node_values_dict[node]
                    
                if i < (len(path_nodes) - 1):
                    next_node = path_nodes[i + 1]
                    if path.get_sign_dict()[(node, next_node)] != "+":
                        node_value = 1.0 - node_value
                        
                    score *= node_value
                else:
                    score *= node_value
                
                # append results to list
                node_values.append(node_value)
                
            # append path results
            path_results.append(score)
        
            print(f'\npath: {path.to_string()}')
            print(f'node values: {node_values}')
            print(f'path score: {path_results}\n')
            print(f'##########################')
        
        # print(f'\n path results: {path_results}')
        result = reduce(lambda x, y: x*y, path_results)

        return f'{hsa}_{elem.get_subpathway().get_start()}_{elem.get_subpathway().get_end()}\t{elem.get_sample()}\t{result}'
