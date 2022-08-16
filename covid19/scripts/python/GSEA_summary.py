#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Ismael Henarejos Castillo
# ihc.europa@gmail.com

import pprint as p

'''
Summary GSEA results by counting functions and if they are up or down
'''


################################ DIRECTORIES ##################################

dir_in = "/home/sinh/workspace/local/covid19/results/6 - modules/"
name = "gsea_covid_summary"


input0 = open(dir_in + name + ".tsv", "r")
print(input0)

################################ OUTPUTS ####################################

output1 = open(dir_in + name + "_done.tsv", "w")

d_fun = {k + "_up" if float(v) >= 0 else k + "_down": 1 for k, v in (((line.replace("\n", "")).split("\t")) for line in input0.readlines())}

input0.close()
input0 = open(dir_in + name + ".tsv", "r")
for line in input0.readlines():
    fun, nes = (line.replace("\n","").split("\t"))
    funup = fun+"_up"
    fundown = fun + "_down"

    if funup in d_fun.keys():
        d_fun[funup] +=  1
    elif fundown in d_fun.keys():
        d_fun[fundown] +=  1

for k,v in d_fun.items():
    output1.write("{}\t{}\n".format(k, v))
