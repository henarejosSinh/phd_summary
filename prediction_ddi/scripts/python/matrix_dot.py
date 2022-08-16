#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Ismael Henarejos Castillo
# ihc.europa@gmail.com

import numpy as np
import time
import pandas as pd
import multiprocessing
from multiprocessing import Process, Array, Queue

# read similarity matrixes
# simmat = pd.read_csv("dot_test_sim_m_all.csv", index_col=0)
# arr = simmat.to_numpy()


def pd_fill_diagonal(df_matrix, value=0):  # fill diagonal with preferred value
    mat = df_matrix.values
    n = mat.shape[0]
    mat[range(n), range(n)] = value
    return pd.DataFrame(mat)

def eval_row_col_product(m1_row, m2_col): 
    '''
    Will perform product of matrices and return max element-wise value 
    '''

    values = [b if a else 0 for a, b in zip(m1_row, m2_col)]
    return(max(values))

def save_max_values(max_val, row_length, control, resList):
    if control != row_length:
        res = str(max_val) + "\t"
        resList.append(res)
        control += 1
        return(resList, control)

    elif control == row_length:  # reached length of row
        res = str(max_val) + "\n"
        resList.append(res)
        control = 0
        return(resList,control)

def matrix_to_dataframe(txt, res_name):

    matrix_dot = np.loadtxt(txt)
    df = pd.DataFrame(
        matrix_dot, index=input1.index.values, columns=input1.columns.values)

    df.to_csv(res_name)

################################### DIRS #######################################

# M1
dir_in = "/home/sinh/workspace/local/ddi/results/03.5_REF_M1/"

# M2
dir_in2 = "/home/sinh/workspace/local/ddi/results/06_vilar_ADE/"
dir_in2 = "/home/sinh/workspace/local/ddi/results/07_kegg/"
dir_in2 = "/home/sinh/workspace/local/ddi/results/04_vilar_chems/"
dir_in2 = "/home/sinh/workspace/local/ddi/results/04_vilar_chems/"
dir_in2 = "/home/sinh/workspace/local/ddi/results/05_vilar_targets,carriers,enzy,transpor/"
dir_in2 = "/home/sinh/workspace/local/ddi/results/03.5_REF_M1/"

# OUT M12

dir_out = "/home/sinh/workspace/local/ddi/results/13_results/"
dir_out = "/home/sinh/workspace/local/ddi/results/09_matrix_product/"

################################### FILES ######################################

### INPUTS

# M1
file1 = dir_in + "M1_art_v2.csv"
file1 = dir_in + "dot_test_m1.csv"  

input1 = pd.read_csv(file1, index_col=0)
arr1 = input1.to_numpy()
print("M1: ", file1)


# M2
file2 = dir_in2 + "ipf_matrix_tc_filled_art.csv"
file2 = dir_in2 + "dot_test_m2.csv"

input2 = pd.read_csv(file2, index_col=0)
arr2 = input2.to_numpy()
print("M2: ", file2)

### Load and print matrixes
M1 = arr1
M2 = arr2
print(M1)
print(M2)

### OUTPUTS 

# MATRIX TXT
res1 = dir_out + "test_base.txt"
output = open(res1, "w")
print("Output (matrix) : ", res1)

# DF CSV
res2 = dir_out + "test_df.csv"
output2 = open(res2, "w")
print("Output (df) : ", res2)

### SEQUENTIAL VERSION

print("starting writing product matrix")
start = time.time()

list_res = []

for i in M1:
    # start with first row
    # the resulting first row of the res matrix will be computed using only the row selected with the rest of the cols
    c = 1 # to count row elements to reconstruct resulting matrix row-wise

    for j in M2.T:  # will iterate over all cols of M2
        max_val = 0
        max_val = eval_row_col_product(i, j)

        # after completing element wise product and retrieving the maximum between all operations of said row and col:
        list_res, c = save_max_values(max_val, len(i), c, list_res)

# Save results in matrix shape
for ele in list_res:
    output.write(ele)
output.close()

# Matrix to df conversion
matrix_to_dataframe(res1, res2)

exit()
