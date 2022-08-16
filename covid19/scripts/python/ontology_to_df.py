#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Ismael Henarejos Castillo 
# ihc.europa@gmail.com

import os

# os.chdir("/home/sinh/Downloads/bingo")
os.chdir("/home/sinh/workspace/local/covid19/results/5 - network_study/bingo")
print(os.getcwd())

# phases "BP", "CC", "MF"
phase = "MF"

# INPUTS
file = "bingo_GOprop5-200_" + phase + ".txt"
file1 = "bingo_gene_GO_" + phase + ".txt"

in0  = open(file, "r")
next(in0)

in1 = open(file1, "r")
next(in1)

# OUTPUTS
out0 = open(phase + "_ont.tsv", "w")
out1= open(phase + "_gene_corr.tsv", "w")

for line in in0.readlines():
    line = line.replace("\n", "")
    corr = line.split(" = ")
    print(corr)
    out0.write("{}\t{}\n".format(corr[0],corr[1]))

for line in in1.readlines():
    line = line.replace("\n", "")
    corr = line.split(" = ")
    print(corr)
    out1.write("{}\t{}\n".format(corr[1],corr[0]))
