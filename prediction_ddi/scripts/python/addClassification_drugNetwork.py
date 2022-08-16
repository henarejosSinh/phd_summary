#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Ismael Henarejos Castillo 
# ihc.europa@gmail.com

from os import write
import pickle
from pprint import *
import pandas as pd
import itertools
import re
from operator import itemgetter

'''
Obtains a SIF and ATT file from DDIs, based on FDR and IVF drugs
'''

################################ DIRECTORIES ####################################

dir = "/home/ihenarejos/workspace/projects/ddi/results/15_manuscript_preparation/"

############################# INPUTS ###########################################

edges = open((dir + "top10IVF_DegreePlus1stN default edge.csv"), "r")

in_class = open("/home/ihenarejos/workspace/projects/ddi/results/15_manuscript_preparation/frd_conditions.tsv", "r") # no header

############################ OUTPUTS ###########################################

outSIF = open(dir + "ivf_sif_2.tsv", "w")

# write headers
outSIF.writelines("drug1\tdrug2\n")

############################## MAIN ############################################

##### 0 Add extra info to the sif

class_frd = {k:v for k,v in ((line.strip("\n")).split("\t") for line in in_class)}
# pprint(class_frd)

check = set()

for line in edges:
    line = line.strip("\n")
    line = line.replace('"', "")
    ddi = line.split(",")[2]
    
    # add relation to the sif
    res = ddi.replace(" (interacts with) ", "\t")
    outSIF.write(res + "\n")
    
    # add each drug with its respective class as a relation
    ddiV2 = ddi.split(" (interacts with) ")
    for i in ddiV2:
        if i in class_frd.keys() and i not in check:
            check.add(i)
            categories = class_frd[i].split(" - ")
                
            for j in categories:
                res = "{}\t{}\n".format(i, j)
                outSIF.write(res)
        
outSIF.close()

