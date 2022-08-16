# Script to detect modules in correlation networks
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Ismael Henarejos Castillo
# ihc.europa@gmail.com

import pprint as p
import re
import os

'''
Filter covid-networks by degree/betweenness/bottlenecks
'''

################################ DEF ####################################

################################ DIRECTORIES ####################################

print(os.getcwd())

# Menstrual cycle phases
os.chdir("/home/sinh/workspace/local/covid19/results/5 - network_study/Mid_Secretory_phase/")
print(os.getcwd())

input0 = "mse_cytohubba_res.csv"
in0 = open(input0, "r")
next(in0)
lines = in0.readlines()

input1 = "mse_network.tsv"
in1 = open(input1, "r")
next(in1)
lines1 = in1.readlines()

################################ OUTPUTS ####################################

output0 = open("mse_top_genes.tsv", "w")
output1 = open("mse_top_sif.tsv", "w")
output1.write("node1\tnode2\tcorr\tpval\tpadj\tabs_corr\n")

################################ MAIN ####################################

########### crete dicts with top genes

# top 20 degree (col 4)
genes = []
top_dict = {}
limit = 9

for index, line in enumerate(sorted(lines, reverse=True, key=lambda line: float(line.split(",")[4]))):
    print(index)
    list_all = line.split(",")
    # list_no_index = "\t".join(list_all)
    res = "{}\t{}".format(list_all[0], "degree")
    genes.append(res)

    if list_all[0] not in top_dict:
        top_dict[list_all[0]] = set()
        top_dict[list_all[0]].add("degree")
    else:
        top_dict[list_all[0]].add("degree")

    if index == limit:
        break


for key,value in top_dict.items():
    string = ",".join(value)
    output0.writelines("{}\t{}\n".format(key, string))

########### READ SIF and create a subSIF with just top genes and COVID-19 genes:
covid = ["ACE2", "TMPRSS2", "TMPRSS4", "MX1", "CTSL", "CTSB", "BSG", "FURIN"]

for ele in covid:
    top_dict[ele] = ""

subsif = [line for line in lines1 if line.replace("\n", "").split(
    "\t")[0] in top_dict and line.replace("\n", "").split("\t")[1] in top_dict]


for i in subsif:
    output1.write(i)
