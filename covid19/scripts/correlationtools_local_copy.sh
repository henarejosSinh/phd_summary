#!/bin/bash
# echo "Tool directory is /data/network/nas01/fivibio-data/projects/JavaProjects/gsrmTools/correlationTool"
#dir "/data/network/nas01/fivibio-data/projects/JavaProjects/gsrmTools/correlationTool"
#echo $dir
# i and o
# echo "Enter absolute path of input file"
#read idir
idir=/data/network/nas01/fivibio-data/projects/covid/results/1-DE/count_matrix_normalised_women.tsv

# echo "Enter absolute path of output file"
# read odir
odir=/data/network/nas01/fivibio-data/projects/covid/results/1-DE/correlation_results.tsv

# echo "Enter number of threads you want to use"
# read $nthreads
nthreads=30

# enter fivibio

# cd $dir
/data/network/nas01/fivibio-data/projects/JavaProjects/gsrmTools/correlationTool/correlationTool correlation -i $idir -o $odir -t $nthreads
