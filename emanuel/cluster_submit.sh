#!/bin/sh
WD=/nfs/nobackup2/saezgrp/homes/emanuel/dream_2014
bsub -n 4 -C 0 -q research-rh6 -M 8000 -R "rusage[mem=8000]" -o $WD/output.txt -e $WD/error.txt "python $WD/$1"
