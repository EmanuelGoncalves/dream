#!/bin/sh
WD=/nfs/nobackup2/saezgrp/homes/emanuel/dream_2014
GENES=$(cat $WD/emanuel/genes.tab)

for GENE in $GENES
do
    bsub -n 4 -C 0 -q research-rh6 -M 12000 -R "rusage[mem=12000]" -o $WD/output.txt -e $WD/error.txt -J $GENE "python $WD/emanuel/leader_sc1_cluster.py $GENE"
done