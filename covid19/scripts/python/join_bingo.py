# Script to detect modules in correlation networks
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Ismael Henarejos Castillo
# ihc.europa@gmail.com

import pprint as p
import re
import os

'''
Join Bingo functional enrichment result tables
'''

################################ DEF ####################################


def check_if_covid(string):
    res = []
    covid_genes = ["ACE2", "TMPRSS2", "TMPRSS4",
                   "MX1", "CTSB", "CTSL", "FURIN", "BSG"]
    if re.search("Genes in test set$", string):
        pass
    else:
        genes = (string.replace("\n", "").split("|"))
        for covid_gene in covid_genes:
            if covid_gene in genes:
                tores = "{}\t".format(covid_gene)
                res.append(tores)
    return(res)


def find_index(lines):
    index = [i for i, item in enumerate(lines) if re.search('GO-ID', item)]
    s = map(str, index)
    index = int(''.join(s))
    no_header = index + 1
    new_list = lines[no_header:]
    return(new_list)


def format_and_return(sublist, output):

    for item in sublist:
        item.replace("\n", "")
        cols = item.split("\t")

        # format GO ID
        max_char = 7
        number_char = len(cols[0]) + 1 # to avoid problem of 0 index
        diff = int(max_char - number_char) + 1
        add = "0" * diff


        cols[0] = "GO:" + add + cols[0]

        list_covid = check_if_covid(cols[-1])  # cols[-1] == function genes
        covid = ",".join([str(s) for s in list_covid])

        # dont add covid for now
        output.write("{}\t{}\t{}\t{}\t{}".format(cols[0],
            cols[1], cols[2], cols[-2], cols[-1]))


################################ DIRECTORIES #################################

print(os.getcwd())

# Menstrual cycle phases
os.chdir("/home/sinh/workspace/local/covid19/results/5 - network_study/Mid_Secretory_phase/")

print(os.getcwd())

################################ OUTPUTS ####################################

output0 = open("mse_processed_join.txt", "w")
output0.write("goid\tpval\tpadj\tdescription\tgenes\n")


################################ MAIN ####################################

r = re.compile("bp|bgo$")

# Add filenames with .bgo extension to a list
filenames = [file for file in os.listdir() if r.search(file)]
filenames.sort()

# For each .bgo file, add formatted results to an unique output file
for file in filenames:
    input0 = open(file, "r")
    lines = input0.readlines()
    sublist = find_index(lines)
    format_and_return(sublist, output0)
