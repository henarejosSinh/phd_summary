#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Ismael Henarejos Castillo 
# ihc.europa@gmail.com

import pickle
import pprint as pp
import pandas as pd
import itertools
import re
from operator import itemgetter

def find_corr(drugs, 
              drug_dict):
    res = set()
    for i in drugs:
        if i in drug_dict:
            res.add(drug_dict[i])

    return res


################################ DIRECTORIES ####################################

dir1 = "/home/ihenarejos/workspace/projects/ddi/results/00_reference_files/"

outdir = "/home/ihenarejos/workspace/projects/ddi/results/14_discussion/"

############################# INPUTS ###########################################

in0 = open(dir1 + "art_192_drugs_14apr21.txt", "r")

in1 = open(outdir + "pca_ddis_revision.txt", "r")

############################ OUTPUTS ###########################################

outname = "interactions_2analyze.tsv"
out = open(outdir + outname, "w")

############################## MAIN ############################################

art_dict = {k:v for k,v in ((line.rstrip("\n").split("\t")[0:2]) for line in in0.readlines())}

lines1 = [i.rstrip("\n") for i in in1.readlines()]

for i in lines1:
    drugs = i.split("_")
    set_corr = find_corr(drugs, art_dict)
    res = " - ".join([s for s in set_corr])

    out.write("{}\n".format(res))

out.close()