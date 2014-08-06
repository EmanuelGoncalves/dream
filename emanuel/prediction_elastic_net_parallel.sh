GENE_LIST=$(cat genes.tab)

# Cluster Paths
CLUSTER_DIR=/nfs/nobackup2/saezgrp/homes/emanuel/dream_2014

for GENE in $GENE_LIST
do
	bsub -n 4 -C 0 -q research-rh6 -M 5000 -R "rusage[mem=5000]" -o $CLUSTER_DIR/output.txt -e $CLUSTER_DIR/error.txt "Rscript prediction_elastic_net_parallel.r $GENE"
done