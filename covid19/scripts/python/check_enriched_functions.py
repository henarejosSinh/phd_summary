#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Ismael Henarejos Castillo
# ihc.europa@gmail.com

import pprint as p
import os
import re

'''
Check number of up/down genes in each enriched function of modules
'''

################################ DEF ###########################################

def check_enriched_function_up(enriched_functions, att_dict):
    function = enriched_functions.split("\t")[1]
    genes = (enriched_functions.split("\t")[-2].split("/"))


    fold_value = 0
    dct = {}
    {k: 1 if float(v) >= fold_value and k not in dct and not dct.update({k:1}) else dct[k] + 1 
    if float(v) >= fold_value and not dct.update({k: dct[k] + 1}) else 1 for (k,v) in ((function,att_dict[gene]) for gene in genes if gene in att_dict.keys() )}

    dct_down = {}
    {k: 1 if float(v) <= fold_value and k not in dct_down and not dct_down.update({k: 1}) else dct_down[k] + 1
     if float(v) <= fold_value and not dct_down.update({k: dct_down[k] + 1}) else 1 for (k, v) in ((function, att_dict[gene]) for gene in genes if gene in att_dict.keys())}

    return(dct, dct_down)

################################ DIRECTORIES ####################################

# dir_in = "/home/sinh/workspace/local/covid19/results/6 - modules/"
dir_in = "/home/ihenarejos/workspace/projects/covid19/results/6 - modules/"
name = "enriched_covid_summary_nofc"

out0 = open(dir_in + name + ".tsv", "w")


# change directory
# os.chdir("/home/sinh/workspace/local/covid19/results/6 - modules")
os.chdir("/home/ihenarejos/workspace/projects/covid19/results/6 - modules")
print(os.getcwd())

# att file
# att = open("/home/sinh/workspace/local/covid19/results/5 - network_study/att_covid.tsv", "r")
att = open("/home/ihenarejos/workspace/projects/covid19/results/5 - network_study/att_covid.tsv", "r")
next(att) # skips header
att = att.readlines()

################################ MAIN ####################################

#### 1 att dict
att_dict = {k:v for (k,v) in ((s.split("\t")[0:2]) for s in att)}
# p.pprint(att_dict)


r = re.compile("ora")

# Add filenames with .bgo extension to a list
filenames = [file for file in os.listdir() if r.search(file)]
filenames.sort()
print(filenames)

out0.write("ID	Name	GeneRatio	BgRatio	pvalue	padjust	qvalue	geneID	Count\tUp\tDown\n")
for file in filenames:
    input0 = open(file, "r")
    next(input0)
    lines = input0.readlines()
    
    for line in lines:
        line = line.replace("\n", "")
        up, down = check_enriched_function_up(line, att_dict)
        res = "{}\t{}\t{}\n".format(line, "".join(
            str(s) for s in [*up.values()]), "".join(str(s) for s in [*down.values()]))
        out0.write(res)

