import numpy as np
import time
import pandas as pd
from contextlib import suppress
import multiprocessing
from multiprocessing import Process, Array, Queue
from contextlib import closing
import ctypes as c

### Sequential version

# t_m2 = pd.read_csv("t_m2_400.txt", sep = "\t")
# print(t_m2.head)
# sim_mat = pd.read_csv("sim_mat_400.txt", sep = "\t")
# print(sim_mat.head)
# start0 = time.time()
# combinations = pd.read_csv("combinations_400.txt", sep = "\t")
# print(combinations.head)
# # print(range(1, len(combinations.columns)))
# end0 = time.time()
# start = time.time()


# shared_arr = Array(c.c_double, len(sim_mat.index) * len(sim_mat.columns))
# arr = np.frombuffer(shared_arr.get_obj())
# arr_res = arr.reshape((len(sim_mat.index), len(sim_mat.columns)))
# print(arr_res)

# # for i in range(0,len(combinations.columns)): # len(combinations.columns)): 
# for i in (combinations.columns):
# 	print(i)
# 	i = int(i.replace("V","")) - 1
# 	print(i)
	
# 	#if  i%100 == 0: # if remainder is zero, print
# 	print("iterator", i)

# 	# print(t_m2.loc[combinations.iloc[0,i],], "\n")
# 	# print(t_m2.loc[combinations.iloc[1,i],])
# 	# criteria1 = t_m2.loc[combinations.iloc[0,i],] != 0
# 	# criteria2 = t_m2.loc[combinations.iloc[0,i],] == t_m2.loc[combinations.iloc[1,i],]
# 	# print(criteria1)
# 	# print(criteria2)
# 	# # print(t_m2.loc[combinations.iloc[0,i],] == 1)
# 	# print(t_m2.loc[combinations.iloc[0,i],] == t_m2.loc[combinations.iloc[1,i],])
# 	nab = sum( (t_m2.loc[combinations.iloc[0,i],] != 0) & (t_m2.loc[combinations.iloc[0,i],] == t_m2.loc[combinations.iloc[1,i],]) )
# 	# criteria_all = criteria1 & criteria2
# 	# print(criteria_all)
# 	# nab = sum(criteria_all)
# 	# method 1 Nab/(Na+Nb+Nab)
# 	# where Na is 1s in A that are not in B and Nb is the inverse
# 	na = sum( (t_m2.loc[combinations.iloc[0,i],] == 1) & (t_m2.loc[combinations.iloc[1,i],] == 0))
# 	nb = sum( (t_m2.loc[combinations.iloc[1,i],] == 1) & (t_m2.loc[combinations.iloc[0,i],] == 0))
# 	# print(nab)
# 	# print(na)
# 	# print(nb)
# 	# print(nab/(na+nb+nab))

# 	# if (na+nb+nab) == 0:  # skip division by zero
# 	#     sim_mat.loc[combinations.iloc[0,i], combinations.iloc[1,i]] = 0
# 	#     sim_mat.loc[combinations.iloc[1,i], combinations.iloc[0,i]] = 0
# 	#     continue
# 	# # 0/x division
# 	# if (nab) == 0:  
# 	#     sim_mat.loc[combinations.iloc[0,i], combinations.iloc[1,i]] = 0
# 	#     sim_mat.loc[combinations.iloc[1,i], combinations.iloc[0,i]] = 0
# 	#     continue

# 	# if everything is okay, fill upper half
# 	# sim_mat.loc[combinations.iloc[0,i], combinations.iloc[1,i]] = (nab/(na + nb + nab))
# 	# try:
# 	with suppress(ZeroDivisionError):

# 	#print(na)
# 	#print(nb)
# 	#print(nab)
# 		# print(nab/(na + nb + nab))
# 		# sim_mat.loc[combinations.iloc[0,i], combinations.iloc[1,i]] = (nab/(na + nb + nab))
# 		print(combinations.iloc[0,i])
# 		print(combinations.iloc[1,i])
# 		print( sim_mat.index.get_loc(combinations.iloc[0,i]))
# 		ri = sim_mat.index.get_loc(combinations.iloc[0,i])
# 		print( sim_mat.columns.get_loc(combinations.iloc[1,i]))
# 		ci = sim_mat.columns.get_loc(combinations.iloc[1,i])
# 		print(type(ri))
# 		print(type(ci))
# 		# save in shared array in corresponding row, col indexes
# 		arr_res[ri,ci] = (nab/(na + nb + nab))

# 		# filling lower half
# 		print( sim_mat.index.get_loc(combinations.iloc[1,i]))
# 		ri2 = sim_mat.index.get_loc(combinations.iloc[1,i])
# 		print( sim_mat.columns.get_loc(combinations.iloc[0,i]))
# 		ci2 = sim_mat.columns.get_loc(combinations.iloc[0,i])
# 		print(type(ri))
# 		print(type(ci))
# 		# save in shared array in corresponding row, col indexes
# 		arr_res[ri2, ci2] = (nab/(na + nb + nab))

# 	# except ZeroDivisionError:
# 	#	pass
# 	# filling the lower half
# 	#sim_mat.loc[combinations.iloc[1,i], combinations.iloc[0,i]] = (nab/(na + nb + nab))
# 	#except:
# 	#	pass

            

# end = time.time()
# print(arr_res)

# df = pd.DataFrame(arr_res, index=sim_mat.index.values, columns=sim_mat.columns.values)
# print(df)
# print("read combinations: ", end0 - start0)
# print("sim_mat filled: ", end - start)

# df.to_csv("sim_mat_filled_py.csv")


### Parallel version

start0 = time.time()

### input files
dir_in = "/home/sinh/workspace/local/ddi/results/00_reference_files/"
dir_in2 = "/home/sinh/workspace/local/ddi/results/07_kegg/"

input1 = dir_in2 + "t_m2_KEGG.txt"  # kegg
input2 = dir_in + "sim_mat_a.txt" 
input3 = dir_in + "combinations_a.txt"

t_m2 = pd.read_csv(input1, sep = "\t")
sim_mat = pd.read_csv(input2, sep = "\t")
combinations = pd.read_csv(input3, sep = "\t")

print(t_m2.head)
print(sim_mat.head)
print(combinations.head)

#### output files

output1 = dir_in2 + "sim_mat_kegg.csv"

end0 = time.time()
start = time.time()

# define functions

def pd_fill_diagonal(df_matrix, value=0):  # fill diagonal with preferred value
    mat = df_matrix.values
    n = mat.shape[0]
    mat[range(n), range(n)] = value
    return pd.DataFrame(mat)

class ProcessAlpha(Process):

	def __init__(self, sim_mat, fingerprints, combinations, batch, array_res):
		Process.__init__(self)
		self.sim_mat = sim_mat
		self.fingerprints = fingerprints
		self.combinations = combinations
		self.batch = batch
		self.array_res = array_res

	def run(self):
	
	# for i in range(0,len(combinations.columns)): # len(combinations.columns)): 
		for i in (self.combinations[self.combinations.columns[self.batch[0]:self.batch[1]]]): # iterate over columns
			# of combinations dataframe (pair of drugs)
			i = int(i.replace("V","")) - 1

			# method 1 Nab/(Na+Nb+Nab)
			# where Na is 1s in A that are not in B and Nb is the inverse
			nab = sum( (self.fingerprints.loc[self.combinations.iloc[0,i],] != 0) & (self.fingerprints.loc[self.combinations.iloc[0,i],] == self.fingerprints.loc[self.combinations.iloc[1,i],]) )
			# translation: where values of A are different than 0 and values of A == B (the order doesn't matter )
			na = sum( (self.fingerprints.loc[self.combinations.iloc[0,i],] == 1) & (self.fingerprints.loc[self.combinations.iloc[1,i],] == 0))
			# where A values are 1 but are 0 in B
			nb = sum( (self.fingerprints.loc[self.combinations.iloc[1,i],] == 1) & (self.fingerprints.loc[self.combinations.iloc[0,i],] == 0))
			# where B values are 1 but are 0 in A

			with suppress(ZeroDivisionError):  # avoid stopping execution in
				# 0/0 cases

				# fill upper half
				# here, retrieve row, col index using combinations
				# combinations.iloc retrieves 
				ri = self.sim_mat.index.get_loc(self.combinations.iloc[0,i])
				ci = self.sim_mat.columns.get_loc(self.combinations.iloc[1,i])
				# save in shared array in corresponding row, col indexes
				self.array_res[ri,ci] = (nab/(na + nb + nab))

				# filling lower half
				ri2 = self.sim_mat.index.get_loc(self.combinations.iloc[1,i])
				ci2 = self.sim_mat.columns.get_loc(self.combinations.iloc[0,i])
				# save in shared array in corresponding row, col indexes
				self.array_res[ri2, ci2] = (nab/(na + nb + nab))


def process_sim_mat(sim_mat, fingerprints, combinations, number_p):

	shared_arr = Array(c.c_double, len(sim_mat.index) * len(sim_mat.columns))
	arr = np.frombuffer(shared_arr.get_obj())
	arr_res = arr.reshape((len(sim_mat.index), len(sim_mat.columns)))
	print(arr_res)

	batch = []
	divs = len(combinations.columns) // number_p

	for i in range(0, number_p):
		batch.append([i * divs, divs * (1+ i)])

	print(batch)
	batch[-1][1] = len(combinations.columns)
	print(batch)
	# print(len(combinations.columns[0:3].tolist()))
	# print(combinations[combinations.columns[0:3]])
	
	# for i in (combinations[combinations.columns[0:3]]):
	# 	if i is not None:
	# 		print(i)
	# 		print(i.replace("V",""))
	# 		print(type(i))

	# start process

	process_list = []
	for i in range(0, number_p):

		print('starting process', i+1)
		process_list.append(ProcessAlpha(sim_mat=sim_mat, 
		fingerprints=fingerprints, combinations=combinations, 
		batch=batch[i], array_res=arr_res))
		process_list[i].start()

	for i in range(0, number_p):
		process_list[i].join()
		print('process ', i + 1, "stopped")

	return(arr_res)


# MAIN
number_p = 6
res = process_sim_mat(sim_mat=sim_mat, fingerprints=t_m2, 
combinations=combinations, number_p=number_p)

end = time.time()
print(res)

df = pd.DataFrame(res, index=sim_mat.index.values, 
columns=sim_mat.columns.values)
df = df.fillna(0)  # nan = 0
print(df)

# diag = 0
df2 = pd_fill_diagonal(df, 0)
df2.index = sim_mat.index.values
df2.columns = sim_mat.columns.values
print(df2)
df2.to_csv(output1)


print("read combinations: ", end0 - start0)
print("sim_mat filled: ", end - start)
