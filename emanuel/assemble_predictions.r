# Set working directory
# wd  <- '~/Documents/projects_data_analysis/dream/' 
wd  <- '/nfs/nobackup2/saezgrp/homes/emanuel/dream_2014/'

setwd(wd)
source(paste(wd, 'dream2014_Import_Datasets.R', sep=''))

# Set-up
genes <- colnames(achilles.training)
samples <- rownames(cn.leader.raw)

gene.files <- list.files(path=paste(wd, 'predictions/', sep=''), full.names=T)

print(gene.files)

predictions <- do.call(rbind, lapply(gene.files, function(file) {
	gene.pred <- read.table(file, sep='\t', header=T, row.names=1)
}))

print(dim(predictions))
predictions <- predictions[genes, samples]
print(dim(predictions))

write.res.gct(predictions, rownames(predictions), colnames(predictions), outfile=paste(wd, 'submissions/umebi_emanuel_21.gct',sep=''))