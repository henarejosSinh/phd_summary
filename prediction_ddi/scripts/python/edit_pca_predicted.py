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

'''
Edit pca results (predicted) to add to paper as suppl
'''

def openSkipheader(file):
    f = open(file, "r")
    next(f)
    return f

def dropFromList(list_in, indexToDrop):
    unwanted = set(indexToDrop)
    list_out = []
    for i in list_in:
        res = [j for j in (i.split("\t")) if i.index(j) not in unwanted]
        res = "\t".join(res)
        list_out.append(res)

    return list_out

################################ Functions ####################################

def name_correspondence(list_in, dict_in0, dict_in1, output_list):

    for line in list_in:

        # ddi modify
        interaction = (line.split("\t")[0]).split("_")
        interactionNames = list(map(lambda x: dict_in0[x], interaction))

        # IVF name
        drugPresence = list(filter(lambda x: x in dict_in1.keys(), interaction))
        drugNamesIVF = list(map(lambda x: dict_in1[x], drugPresence))

        line = line.split("\t")
        line[0] = "-".join(interactionNames)
        line[5] = "-".join(drugNamesIVF)
        
        output_list.append(("\t".join(line)))
        
    return output_list

################################ INPUTS ####################################

# PCA RESULT FILES
pcaPred = "/home/ihenarejos/workspace/projects/ddi/results/15_manuscript_preparation/ddis_results - PCA_V2.tsv"

# DRUGBANK ID-DRUG NAME CORRESPONDENCE

DBcorrIvf = "/home/ihenarejos/workspace/projects/ddi/results/00_reference_files/art_192_drugs_14apr21.txt"
DBcorrIvf = {k:v for (k,v) in ((line.strip("\n")).split("\t") for line in open(DBcorrIvf, "r"))}

DBcorr = "/home/ihenarejos/workspace/projects/ddi/results/00_reference_files/drugbank_info_2020_07.txt"
DBcorr = {k:v for (k,v) in ((line.strip("\n")).split("\t")[0:2] for line in open(DBcorr, "r"))}

################################ OUTPUT ####################################

pcaPredFiltered = "/home/ihenarejos/workspace/projects/ddi/results/15_manuscript_preparation/pca_pred.tsv"

################################ MAIN ####################################

pred = [line.strip("\n") for line in openSkipheader(pcaPred)]

# indexDrop = [2,3,4,6,7,8,9,14]
# pred = dropFromList(pred,indexDrop)
# print(pred[1:2])
# exit()

output_list = []
pred_filtered = (name_correspondence(pred, DBcorr, DBcorrIvf, output_list))

# filter columns and fix associated effect name change

output = open(pcaPredFiltered, "w") 
output.write("ddi\tscore\tfrd\tpredicted effect\tatc codes\tguidelines\tclinical trials\tinteraction type\n")

unwanted = {2,3,4,6,7,8,9,14}
for i in pred_filtered:

    # format old columns to show name of drug
    oldColumns = i.split("\t")
    d1 = oldColumns[0].split("-")[0]
    d2 = oldColumns[0].split("-")[1]
    
    # atc
    atc = oldColumns[11].split("_")
    if atc[0] == "NA":
       oldColumns[11] = "{}:({})-{}:({})".format(d1, "NA", d2, "NA")
    else:
        oldColumns[11] = "{}:({})-{}:({})".format(d1, atc[0], d2, atc[1])

    # guidelines
    guide = oldColumns[12].split("_")
    oldColumns[12] = "{}:({})-{}:({})".format(d1, guide[0], d2, guide[1])

    # clinical trials
    trials = oldColumns[13].split("_")
    oldColumns[13] = "{}:({})-{}:({})".format(d1, trials[0], d2, trials[1])

    # fix predicted effect names
    predEffect = oldColumns[10]
    oldDrug1, oldDrug2 = oldColumns[8].split("_")

    if re.search(d1, predEffect):
        oldColumns[10] = re.sub(DBcorr[oldDrug2], d2, predEffect)
    if re.search(d2, predEffect):
        oldColumns[10] = re.sub(DBcorr[oldDrug1], d1, predEffect)

    # filter unwanted columns
    newList = [i for i in oldColumns if oldColumns.index(i) not in unwanted]

    output.write("\t".join(newList) + "\n")
