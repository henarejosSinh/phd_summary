# ihc.europa@gmail.com
# Functions for pva KEGG project

import re
import traceback
from collections import namedtuple

# classes from project
from AttEntry import *

# ========================
# functions
# ========================

def parseEffectScoreFile(effect_file):
    '''
    # parses the variant effect score file and saves info to a set
    returns : set refScores
    '''
    refScores = set()
    
    try:
        lineList = [line.strip("\n") for line in effect_file]
        ref_score = namedtuple( 'reference_score', ["effect", "impact", "score"])
        for i in lineList:
            effect, impact, score = i.split("\t")
            refScore = ref_score(effect, impact, score)
            refScores.add(refScore)
            
    except Exception:
        traceback.print_exc()
        print(f'Error while parsing variant effect scores file: {effect_file}')
        
    return refScores

def parseAttributeFile(att_dir, att_file):
    '''
    parses kegg attribute file and returns a dictionary with 
    att_entry_class objects
    Does not save xcord and ycord for each entry
    returns: dict att_entries
    '''
    
    try:
        att_entries = {}
        file = open(att_dir+att_file, 'r')
        next(file) # skips header
        lineList = [line.strip("\n") for line in file]
        for line in lineList:
            nodeID, genes, types, genesNames, nodeName, *_ = line.split("\t")
            genes = genes.replace("\|", "||").replace("\&", "&&")
            genesNames = genesNames.replace("\|", "||").replace("\&", "&&")
            att_entries[nodeID] = AttEntry.create(nodeID, genes, types, genesNames,nodeName)
            
    except Exception as e:
        print(f'\nERROR parsing attribute file: {att_file} ')
        traceback.print_exc()
        
    return att_entries

def evaluator(text: str) -> float:
    sp1 = (text.split(' & '))

    d = {}
    for i in range(0, len(sp1)):
        elem = re.sub('[\(\)]', '', sp1[i])
        if ' | ' in elem:
            value = max(elem.split(' | '))
            d[i] = value
        else:
            d[i] = elem
    
    # print(text.split(' & '))
    # print(d)
    # print(min(d.values()))
    
    return float(min(d.values()))

# ========================
# test code here
# ========================

# # boolean eval 
# import re
# test = '(ACTC1) & (TPM4 | TPM3 | TPM1 | TPM2)'
# test = '(0.8) & (0.2 | 0 | 0.1 | 0.5)'
# test = '(TAS1R2) & (TAS1R3)'
# test = '(0.1) & (0.2)'
# test = 'ADCY9 | ADCY6 | ADCY3'
# test = '0.0000005 | 0.2 | 1'

# # str = "[Hi all], [this is] [an example] "
# contents = re.findall('\((.*?)\)', test)
# # res = (0.8) * (max(0.2, 0, 0.1, 0.5))
# # print(res)
# # print(test.split('('))
# print(contents)
# for i  in contents:
#     print(i)

# # & 
# print(test.split(' & '))
# sp1 = (test.split(' & '))

# d = {}
# for i in range(0, len(sp1)):
#     elem = re.sub('[\(\)]', '', sp1[i])
#     if ' | ' in elem:
#         value = max(elem.split(' | '))
#         d[i] = value
#     else:
#         d[i] = elem
# pprint(d)
# print(min(d.values()))

# # orlist = []
# # for i in sp1:
# #     if '|' in i:
# #         print(i.split(' | '))
# #         res = i.split(' | ')
# #         orlist.append(res)
        
# # print(orlist)

# # for i in orlist:
# #     min_or = max(i)
    
# # print(min_or)