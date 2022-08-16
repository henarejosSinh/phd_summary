#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Ismael Henarejos Castillo 
# ihc.europa@gmail.com

import numpy as np
import time
import pandas as pd
from contextlib import suppress
import multiprocessing
from multiprocessing import Process, Array, Queue
from contextlib import closing
import ctypes as c

# transform results from openbabel/ipf|ade/targets vectors to matrix format

start = time.time()

# files to use
# file_name = "/home/sinh/workspace/local/ddi/results/06_vilar_ADE/tc_ade.csv"  #ade
file_name = "/home/sinh/Downloads/08_interactome/interactome_tc_v2.csv"  # interactome

tc_file = open(file_name, 'r')
next(tc_file) # skips header of tc file (drug1,drug2,tc)
lines = tc_file.readlines()
# sim_mat = pd.read_csv("/home/sinh/workspace/local/ddi/results/00_reference_files/sim_mat_a.txt", sep = "\t")
sim_mat = pd.read_csv("/home/sinh/Downloads/sim_mat_a.txt", sep = "\t")

# out
res_dir = "/home/sinh/workspace/local/ddi/results/03.5_REF_M1/"
res_dir = "/home/sinh/Downloads/08_interactome/"

res_file = res_dir + "interactome_matrix_tc_filled_art.csv"

# define functions
def pd_fill_diagonal(df_matrix, value=0): #fill diagonal with preferred value
    mat = df_matrix.values
    n = mat.shape[0]
    mat[range(n), range(n)] = value
    return pd.DataFrame(mat)

# define class

class Processbeta(Process):

    def __init__(self, sim_mat, tc_file, batch, array_res):
        Process.__init__(self)
        self.sim_mat = sim_mat
        self.tc_file = tc_file
        self.batch = batch
        self.array_res = array_res

    def run(self):

    # c = 1 # index of lines
        for line in self.tc_file[self.batch[0]:self.batch[1]]:
            # print(c)
            line = line.replace("\n", "")
            ddi = line.split(",")

            # GET D1, D2, TC values
            d1 = ddi[0]
            d2 = ddi[1]
            tc = ddi[2]

            # use pd.loc to search in rownames/colnames of the matrix
            # fill upper/lower half of matrix correspondingly 
            if d1 in self.sim_mat.index.values \
                and d2 in self.sim_mat.index.values:

                # here, retrieve row, col index in sim_mat
                ri = self.sim_mat.index.get_loc(d1)
                ci = self.sim_mat.columns.get_loc(d2)
                # save in shared array in corresponding row, col indexes
                self.array_res[ri,ci] = tc

                # filling lower half
                ri2 = self.sim_mat.index.get_loc(d2)
                ci2 = self.sim_mat.columns.get_loc(d1)
                # save in shared array in corresponding row, col indexes
                self.array_res[ri2, ci2] = tc

                # c = c + 1

def process_sim_mat(sim_mat, tc_file, number_p):

    shared_arr = Array(c.c_double, len(sim_mat.index) * len(sim_mat.columns))
    arr = np.frombuffer(shared_arr.get_obj()) # allow use of numpy
    arr_res = arr.reshape((len(sim_mat.index), len(sim_mat.columns)))
    print(arr_res) # where tc values will bet stored, having the same 
    # dim as sim_mat

    batch = []
    num_lines = sum(1 for line in tc_file)
    print(num_lines)
    divs = num_lines // number_p

    for i in range(0, number_p):
        batch.append([i * divs, divs * (1+ i)])

    print(batch)
    batch[-1][1] = num_lines # ensures last process get too work until last line
    print(batch)

    # print(lines[batch[0][0]:batch[0][1]])

    # start process

    process_list = []
    for i in range(0, number_p):

        print('starting process', i+1)
        process_list.append(Processbeta(sim_mat=sim_mat, tc_file=lines, 
        batch=batch[i], array_res=arr_res))
        process_list[i].start()

    for i in range(0, number_p):
        process_list[i].join()
        print('process ', i + 1, "stopped")

    return(arr_res)

# MAIN
number_p = 6
res = process_sim_mat(sim_mat=sim_mat, tc_file=lines, number_p=number_p)

end = time.time()
print(res)

# create dataframe from matrix using labels (values) of rownames/colnames
df = pd.DataFrame(res, index=sim_mat.index.values, columns=sim_mat.columns.values)
df = df.fillna(0) # nan = 0
print(df)

# diag = 0
df2 = pd_fill_diagonal(df, 0)
df2.index = sim_mat.index.values
df2.columns =sim_mat.columns.values
print(df2)
df2.to_csv(res_file)

end = time.time()
timer = end-start

print("matrix conversion %s s \n" %(timer))
