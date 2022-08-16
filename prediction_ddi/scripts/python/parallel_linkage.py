#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Ismael Henarejos Castillo
# ihc.europa@gmail.com

import numpy as np
import time
import pandas as pd
import multiprocessing
from multiprocessing import Process, Array, Queue

'''
Script to obtain a list of linked ddis based on a tc cut off value.
It will retrieve the ddi, the max value obtained through matrix product and
the tc value
'''


def pd_fill_diagonal(df_matrix, value=0):  # fill diagonal with preferred value
    mat = df_matrix.values
    n = mat.shape[0]
    mat[range(n), range(n)] = value
    return pd.DataFrame(mat)


def eval_row_col_product_retrieve_index(m1_row, m2_col):
    '''
    Will perform product of matrices and return max element-wise value 
    '''

    values = [b if a else 0 for a, b in zip(m1_row, m2_col)]
    max_val = (max(values))
    index = m2_col.tolist().index(max_val)
    # Note: The index() method only returns the first occurrence of the matching element.
    return(max_val, index)

### REDEFINE. IS FAULTY
# i.e 2675,3345, (two index at some point)
def save_linkage(max_val, label, list_labels, row_index, col_index, control, resList, batch_end, row_length):

    if control != row_length:

        if max_val == 0 and col_index != row_index:
            res = "{}_{},{},{}\n".format(
                list_labels[row_index], list_labels[col_index], "", 0)
            control = control + 1
            resList.append(res)
            if row_index != batch_end:
                    resList.append("{},".format(row_index))

        else:
            max_rel = '{}_{}'.format(
                list_labels[row_index], label)

            if max_val >= 0 and col_index != row_index:

                res = "{}_{},{},{}\n".format(
                    list_labels[row_index], list_labels[col_index], max_rel, max_val)
                control = control + 1
                resList.append(res)

                if row_index != batch_end:
                    resList.append("{},".format(row_index))
    
    elif control == row_length:

        if max_val == 0 and col_index != row_index:
            res = "{}_{},{},{}\n".format(
                list_labels[row_index], list_labels[col_index], "", 0)
            control = 0
            resList.append(res)
        if row_index != batch_end:
            resList.append("{},".format(row_index))

        else:
            max_rel = '{}_{}'.format(
                list_labels[row_index], label)

            if max_val >= 0 and col_index != row_index:

                res = "{}_{},{},{}\n".format(
                    list_labels[row_index], list_labels[col_index], max_rel, max_val)
                control = 0
                resList.append(res)
                if row_index != batch_end:
                    resList.append("{},".format(row_index))

    return(resList, control)


def matrix_to_dataframe(txt, res_name):

    matrix_dot = np.loadtxt(txt)
    df = pd.DataFrame(
        matrix_dot, index=input1.index.values, columns=input1.columns.values)

    df.to_csv(res_name)

### Parallel Version

# modify process class


class ProcessOmega(Process):

    def __init__(self, matrix1, matrix2, batch, cola, labels, index):
        Process.__init__(self)
        self.cola = cola
        self.matrix1 = matrix1
        self.matrix2 = matrix2
        self.batch = batch
        self.labels = labels
        self.index = index

    def run(self):

        if self.batch[0] == self.batch[1]:
            exit()
        
        row_index = self.batch[0]
        list_res = []
        list_res.append("{},".format(row_index))

        for i in self.matrix1[self.batch[0]:self.batch[1], ]:
            # start with first row
            # the resulting first row of the res matrix will be computed using only the row selected with the rest of the cols
            # print(i)
            c = 1  # to count row elements to reconstruct resulting matrix row-wise
            x = 1  # to control col_index

            for j in self.matrix2.T:  # will iterate over all cols of M2

                max_val = 0
                max_val, index = eval_row_col_product_retrieve_index(i, j)
                
                label = self.labels[index]
                col_index = x - 1

                list_res, c = save_linkage(
                    max_val, label,  self.labels, row_index, col_index, c, list_res, self.batch[1], len(i))
                # col_index control
                x = x + 1
            # next row index
            row_index = row_index + 1
            

        self.cola.put(list_res)

# RESTET MAX REL SPECIALLY COLS


def dot_process(M1, M2, num_process, labels):

    cola = multiprocessing.Manager().Queue()
    divs = (M1.shape[0]) // num_process  # number of rows divided by process
    batch = []

    print(M1.shape[0])
    print(divs)

    for i in range(0, num_process):
        batch.append([i*divs, (i + 1) * divs])

    batch[-1][1] = M1.shape[0]
    print(batch)
    process_list = []

    for i in range(0, num_process):

        print('Launching process: ' + str(i + 1))

        process_list.append(ProcessOmega(
            M1, M2, batch[i], cola, labels, i))
        process_list[i].start()

        # while not cola.empty():
        #     lists_of_res=cola.get()
        #     print(lists_of_res)
        #     res.append(lists_of_res)

    c = 0
    for i in process_list:
        c += 1
        i.join()
        print('Process number ', c, ' stopped.')

    res = []
    while not cola.empty():
        lists_of_res = cola.get()
        res.append(lists_of_res)

    return(res)


################################ DIRECTORIES ####################################

root_dir = "/home/sinh/workspace/local/ddi/results/"

# # M1
dir_in = root_dir + "03.5_REF_M1/"

# # M2
m2folder = "08_interactome/"
dir_in2 = root_dir + m2folder


# Linkage results
dir_out = root_dir + "13_results/"

############################# INPUTS ###########################################

# M1
file1 = dir_in + "dot_test_m1.csv"
file1 = dir_in + "M1_art_v2.csv"

# M2
csv = "dot_test_m2.csv"
csv = "interactome_matrix_tc_filled_art.csv"

file2 = dir_in2 + csv

# OPEN FILES

# M1
input1 = pd.read_csv(file1, index_col=0)
arr1 = input1.to_numpy()
print("M1: ", file1)

# M2
input2 = pd.read_csv(file2, index_col=0)
arr2 = input2.to_numpy()
print("M2: ", file2)

### PRINT MATRIX
M1 = arr1
M2 = arr2
print(M1)
print(M2)

############################ OUTPUTS ###########################################

txt = "int_linkage.csv"
res1 = dir_out + txt

output = open(res1, "w")
output.write("index,ddi,max_interaction,maxVal\n")
print("Output (list) : ", res1)

############################## MAIN ############################################

start = time.time()
num_process = 6
labels = input1.index.values  # extract labels to use in the linkage

res = dot_process(M1, M2, num_process, labels)

list_res = [item for sublist in res for item in sublist]

for ele in list_res:
    output.write(ele)
output.close()

end = time.time()
print("linkage list written in: {} s ".format(end - start))
