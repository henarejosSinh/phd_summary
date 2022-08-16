# ihc.europa@gmail.com
# count drugs in tsv files if involved in certain interactions

from pprint import *
file1 = open('/home/ihenarejos/workspace/projects/ddi/Supplemental table 3 - 3C.tsv', 'r')
output = open('/home/ihenarejos/workspace/projects/ddi/res3c', 'w')

counts = {}

for line in file1:
    d1, d2, *_ = line.strip().split("\t")
    if d1 in counts.keys():
        counts[d1] += 1
    else:
        counts[d1] = 1
    
    if d2 in counts.keys():
        counts[d2] += 1
    else:
        counts[d2] = 1
        
res = sorted(counts.items(), key=lambda x:x[1], reverse=True )

for tup in res:
    output.writelines(f'{tup[0]}\t{tup[1]}\n')