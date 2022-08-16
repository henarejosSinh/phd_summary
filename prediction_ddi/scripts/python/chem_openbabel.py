#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import subprocess
import re
import os
import time


# time start

#  2468 drugs x 2467 TC Calculation (aprox 6,091,024 combinations)

start = time.time()

# 2. Specify the TC cutoff. This option is useful if only the TCs of similar molecules above the established cutoff are needed. Otherwise,
# set T_CUTOFF=0 to provide all TC pair values:


T_CUTOFF = 0

# Specify fingerprint (‘FP2’, ‘FP3’, ‘FP4’ or ‘MACCS’); as an example, to specify fingerprint MACCS, use the following command:

FINGERPRINT = 'MACCS'

# Open smiles

res_dir = "/home/sinh/workspace/local/ddi/results/04_vilar_chems/"

filename = res_dir + "smiles_approv_drugbank_noblanks_id.tsv"
filename = res_dir + "smiles_approv_drugbank_noblanks_id_2.tsv"
# filename = res_dir + "smiles_approv_drugbank_test400.csv"
in_file = open(filename, "r")

temp_name = "temp_smi_file.txt"  ### MUST BE IN SAME DIRECTORY AS SCRIPT OR ELSE WONT WORK
input_temp = open(temp_name, 'w')

# dict of drugs to compare

input_dict = {}

for line in in_file:
    newline = (line.replace("\n", "").split('\t'))
    id_name = newline[0]
    smiles = newline[1]

    if smiles != "":

        input_dict[id_name] = smiles
        input_temp.write('%s\t%s\n' %(smiles, id_name) )

in_file.close()
input_temp.close()

# check dict
# print(input_dict)

# open results file
res_file = res_dir + "TC_results.csv"
res_file = res_dir + "TC_results_2.csv"
f = open(res_file, 'w')

writer = csv.writer(f)
writer.writerow(['chemical1', 'chemical2', 'TC'])

# For each chemical in input list, calculate the TC between that chemical and all other chemicals in the input list using Open Babel:

for chemical1 in input_dict:
    # if chemical1 != "":

    babel_command = 'obabel -ismi -:"%s" temp_smi_file.txt -ofpt -xf%s' % (input_dict[chemical1], FINGERPRINT)
    # babel_command = 'obabel -ismi -:"[N]=O" temp_smi_file.txt -ofpt -xfMACCS'
    # print(babel_command)
    output = subprocess.Popen(babel_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')

# Read and parse output from Open Babel:
    # for line in output:
    #     print(line)

# for line in output.stderr:
#     print(line)

    TC_list = []

    for line in output.stdout:
        # print(line)

        if line != '':
        
            newline = re.split('>|=', line)
            #newline: ['', 'CHEMBL1382 Tanimoto from CHEMBL973 ', ' 0.2\n']
            #indices: [0] [1] [2]
            if len(newline) > 2:
                id_catcher = newline[1].split()
                chemical2 = id_catcher[0]
                TC = float(newline[2].strip())
                TC_list.append((chemical2, TC))
        else:
            break

    # print(TC_list)
# 11. Write the TCs exceeding the cutoff to the output file (exclude chemical1 = chemical2—exclude chemicals with the same molecule
#name—where TC = 1):
    for chemical2,TC in TC_list:
        if TC > T_CUTOFF and chemical1 != chemical2:
            writer.writerow([chemical1, chemical2, TC])

f.close()

# time end
end = time.time() # 22jan21 > 2764s = ~46 min

print("tanimoto coeff. calculated in : ", end - start)
# os.remove('temp_smi_file.txt')

# 12. After the script is run in Python, the output file will produce a list of pairs of compounds and the relevant TC that quantifies the
# level of similarity between them. Transform the data thus obtained into a matrix M2 containing a TC in each cell and set values in the
# diagonal to 0, as detailed in Step 10 of the main PROCEDURE.


