# ihc.europa@gmail.com
# script to calculate the score of a given node

# imports
from pprint import pprint
import statistics
import re
from functools import reduce
from abc import ABC, abstractclassmethod, abstractmethod

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
class scoreAnalysis:
    
    def __init__(self, hsa_genes_dict,
                 variants_x_genes, effect_scores: EffectScore, vcf : Vcf, samples):
        
        self.hsa_genes_dict = hsa_genes_dict
        self.variants_x_genes = variants_x_genes
        self.effect_scores = effect_scores
        self.vcf = vcf
        self.population = samples
        
        self.gene_pattern = re.compile('((hsa:)?([a-zA-Z\-0-9]+))')
        self.eval = evaluator
        
    
    def __calculate_score(self, changes):
        '''
        '''
        
        values = []
        
        for effect_score in self.effect_scores:
            if getattr(effect_score, 'effect') in changes:
                #print(f'''
                    #effect score namedtuple: {effect_score}
                    #effect {getattr(effect_score, 'effect')}
                    #score {getattr(effect_score, 'score')}
                    #''')
                values.append(float(getattr(effect_score, 'score')))

        return min(values)  # if there's more than one change possible affecting the same gene in a variant, select the worst effect
            
    def __get_gene_value(self, gene, sample, variants_x_gene): # 2d dict
        
        if gene not in variants_x_gene.keys():
            return 1.0
        gene_score = 1.0 
        
        variants_indexes = variants_x_gene[gene]
        for index in variants_indexes:
            # # print(f'\n{variant.variant}')
                                           
            changes: list
            allele_score: list
            
            allele_score = []
            
            # to avoid RAM problems 
            variant = Variant(self.vcf.list_variants[index], sample_list=self.population)
            
            # set up
            alleles_and_genes = variant.get_snpeff_allele_genes()
            allele_corr = variant.allele_correspondence()

            variant_sample = VariantSample(variant, sample)
            
            # extract the haplotype of a given variant for a certain sample,
            # then retrieve its alleles and check if they are reference or not
            alleles = re.split("\/|\|", variant_sample.haplotype)
            #print(alleles)
            for allele in alleles:
                
                if allele == '0' or allele == '.':
                    allele_score.append(1.0) 
                    continue
                    
                elif allele_corr[allele] in alleles_and_genes.keys():
                    nuc = allele_corr[allele] 
                    annotations = variant.snpEff[nuc]
                    
                    # iterate over transcripts corresponding to alleles(nuc)
                    # and retrieve changes
                    for gene_snpeff in alleles_and_genes[nuc]:
                        if gene_snpeff == gene:
                            changes = [transcrit['Annotation'].split('&') for transcrit in annotations]
                    
                    changes = [item for sublist in changes for item in sublist]
                    
                    # calculate score considering changes -- add all changes
                    score = self.__calculate_score(changes)
                    allele_score.append(score)
                    
            # # print(f'{allele_score}\t{variant_sample.to_string()}')
    
            # calculate mean between two alleles and multiply for the result for each variant affecting the gene
            if len(allele_score) == 0:
                allele_score.append(1.0)
            gene_score *= statistics.mean(allele_score) 
            
        # print(f'\ngene score: {gene_score}') 
        return gene_score
    
    def __calculate_node_value(self, node, att_dict,
                                sample, variants_x_gene):
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
        
        # # print(f'''
        #      .att entry: {attributes.to_string()}
        #       genes in entry: {formula}
        #       matcher: {matcher}
        #       genes: {genes}
        #       ''')
        
        for gene in genes:
            gene_value = self.__get_gene_value(gene, sample, variants_x_gene)
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
                                                             sample, self.variants_x_genes)
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
        
            # # print(f'path: {path.to_string()}')
            # # print(f'node values: {node_values}')
            # # print(f'path score: {path_results}\n')
            # # print(f'##########################')
        
        # # print(f'\n path results: {path_results}')
        result = reduce(lambda x, y: x*y, path_results)

        return f'{hsa}_{elem.get_subpathway().get_start()}_{elem.get_subpathway().get_end()}\t{elem.get_sample()}\t{result}'
    
                    

    

    
    
