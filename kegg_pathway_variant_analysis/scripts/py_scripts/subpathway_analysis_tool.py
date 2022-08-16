# ihc.europa@gmail.com
# main tool of subpathway variant analysis

# IMPORT a bunch of things
from EffectScore import EffectScore
from subpathway_analysis_io_reader import SubpathwayAnalysisInput, SubpathwayAnalysisReader
import sys
import os
import re
import argparse
import csv
from pprint import *
from collections import namedtuple  # to name tuples
from collections import defaultdict # for multidimensional dicts
import traceback
import random

# imports from other scripts
from AttEntry import *
import Path
from Subpathway import *
from functions import *

# import classes from other directories
sys.path.insert(1, '/home/ihenarejos/workspace/projects/pof/scripts/python')
from class_vcf import *
from class_variant import *

# ========================
# Main
# ========================

# pars = argparse.ArgumentParser(description= 'Pva: Pathway variant analysis')
# pars.add_argument("-att", "--attFile-directory", 
#                 help="Directory with all the attribute files",
#                 required=True)
# pars.add_argument("-vcf", "--vcf-file", help='VCF file', required=True)
# pars.add_argument("-eff", "--effect-score-file", 
#                   help='effect score (3-col tab separated file: effect, impact, \
#                       score, with header', required=True)
# pars.add_argument("-hsa", "--hsa-gene-correspondence", help="hsa-gene correspondence",
#                   required=True)
# args = pars.parse_args()

# open files
# if args.vcf is not None:
#     vcf = open(args.vcf, "r")
# if args.eff is not None:
#     eff = open(args.eff, "r")
# if args.att is not None:
#     att = open(args.att, "r")  MUST BE DIR
# if args.hsa is not None:
#     hsa = open(args.hsa, "r") 
#     next(hsa) # since it has header, skip it:

################################################################################


################################ test ##########################################

#vcf = '/home/ihenarejos/workspace/projects/pathway_variant_analysis/results/06_pva_adaption/test_files/hsa04614.vcf'
vcf = '/home/ihenarejos/workspace/projects/pathway_variant_analysis/results/03_pva_score/snpSift_1000g_annot_mod_to_missings_qual.vcf'
# vcf = open('/home/ihenarejos/workspace/projects/pathway_variant_analysis/results/06_pva_adaption/snpSift_1000g_annot_missing_sift.vcf', 'r')

eff = open('/home/ihenarejos/workspace/projects/pathway_variant_analysis/results/06_pva_adaption/test_files/bio_score.txt', 'r')
next(eff)

#att = '/home/ihenarejos/workspace/projects/pathway_variant_analysis/results/06_pva_adaption/test_files/att/'
att = '/home/ihenarejos/workspace/projects/pathway_variant_analysis/results/06_pva_adaption/ATTs/'

hsa = open('/home/ihenarejos/workspace/projects/pathway_variant_analysis/results/06_pva_adaption/test_files/corresponde_vcf_hsa_with_syns.tsv', 'r')

################################################################################

# dir = os.getcwd()
#dir = '/home/ihenarejos/workspace/projects/pathway_variant_analysis/results/06_pva_adaption/test_files/subpathway_tests/'
dir = '/home/ihenarejos/workspace/projects/pathway_variant_analysis/results/06_pva_adaption/results/'
outdir = '/home/ihenarejos/workspace/projects/pathway_variant_analysis/results/06_pva_adaption/test_files/outdir'
# outfile_name = '/out_test.txt'
outfile_name = '/out_pva_kegg.txt'


# run, initialize variables (reminder of variables and type to be)
subPathDir = dir
subPathFiles = []
outputPath = outdir
variantList = []
hsaGeneDict = {}  # string, string (hsa, gene)
sampleNames = set() 
hsaSubpaths = {} # hsa, subpaths
variants_x_gene = {} # gene, variant
attPathFiles = []
attFileDict = {} # id, path

# compile patter to get hsa
pattern = re.compile("(hsa\\d{5})")

# read hsa correspondence file and fill up dictionary with it
hsaGeneDict = {"hsa:{}".format(hsaId):gene for gene, hsaId in 
               ((i.strip("\n").split("\t"))[0:2] for i in hsa)}

# prepare variant object from vcf file
vcfFile = Vcf(vcf)

# =========================================
# add variant list to dict by gene affected
# =========================================
# print(f'\nAdding variants to dict')

if len(vcfFile.list_variants) > 0:
    #print(f'number of variants:')
    vcfFile.num_variants()
    sampleNames = vcfFile.sample_list
    # iterate and get genes from variants affecting protein coding genes:
    for i in range(0, len(vcfFile.list_variants)):
    #for i in range(0, 10):
        genes = Variant(vcfFile.list_variants[i], sample_list=sampleNames).get_genes_in_variant() # unique list
        for j in genes:
            if j not in variants_x_gene.keys(): 
                variants_x_gene[j] = []
                variants_x_gene[j].append(i)
                # variants_x_gene[j].append(Variant(i, sample_list=sampleNames))
            else:
                variants_x_gene[j].append(i)
                # variants_x_gene[j].append(Variant(i, sample_list=sampleNames))
                

# print(f'Variants added to dict: {list(variants_x_gene.items())[0:2]}')

# ======================
# creating list of files
# ======================

#print(f'\nReading subpathway and att files')

try:
    for *_, filenames in os.walk(subPathDir):
        subPathFiles.extend(filenames)
except:
    print(f'Error when loading and reading subpathway files')
    traceback.print_exc()

try: 
    for *_, filenames in os.walk(os.path.abspath(att)):
        attPathFiles.extend(filenames)
        
        for file in attPathFiles:
            # add files to att dict with hsa as id
            matcher = pattern.match(file)
            if matcher:
                attFileDict[matcher.group(1)] = file
             
except:
    print(f'Error when loading and reading .att kegg files')
    traceback.print_exc()

# print(f'Subpathway files: {subPathFiles[0:5]} \t .att files: {list(attFileDict.items())[0:5]}\n')

# ========================
# parse effect, att files*
# ========================
# print(f'\nParsing effect and att files')

refScores = parseEffectScoreFile(eff)
# print(f'\nRef scores file {random.sample(refScores, 5)}\n')
# # pprint(refScores)

# parse attribute using subpath files
for file in subPathFiles:

    matcher = pattern.match(file)
    if matcher:
        hsa = matcher.group(1)
        if hsa not in attFileDict.keys():
            print(f'.ATT file not found for kegg id:\n{hsa}')
            exit()
        
        # if successful, parse corresponding att file:
        attEntriesMap = parseAttributeFile(att, attFileDict[hsa])
        
        # read corresponding subpath file
        read_file = open(dir+file, 'r')
        next(read_file)
        
        dict_2d = defaultdict(dict)
        try:
            for row in read_file:
                row = row.strip("\n")
                # print(f'\n{row}\n')
                
                start, end, _, path_se = row.split("\t")
                p = Path.Path(start, end, path_se)
                if not dict_2d.get(start, {}).get(end, {}):
                    subpath = Subpathway(start, end, hsa, attEntriesMap)
                    dict_2d[start][end] = subpath
                else:
                    subpath = dict_2d.get(start, {}).get(end, {})
                    
                # add path to subpathway 
                subpath.add_path(p)
            
            
            # put information in the hsa, subpaths dict:
            # # pprint(dict_2d)
            # for sp in dict_2d.values():
                # hsaSubpaths[hsa] = sp
            hsaSubpaths[hsa] = dict_2d
            
        except:
            print(f'\nERROR; failed at creating multidimensional dictionary')
            traceback.print_exc()
            print(file)
            exit()

# # print(f'hsa_subpath dict: {list(hsaSubpaths.items())[0:5]} \t att dict: {list(attEntriesMap.items())[0:5]}\n')
# pprint(f'\n{dict_2d}')
# print(f'\nhsa_subpath dict: {(hsaSubpaths)} \t att dict: {list(attEntriesMap.items())}')
      
        
##########################################################
### IMPLEMENT iterative or parallel execution of analysis
##########################################################

output = open(outputPath + outfile_name, 'w')
output.writelines(f'id\tsample\tresult')

reader = SubpathwayAnalysisReader(hsaSubpaths, sampleNames, variants_x_gene, refScores, vcfFile)
# print(reader.hsa_subpathways)

out = reader.run()

for i in out:
    res = i.get_result()
    output.writelines(f'\n{res}')
