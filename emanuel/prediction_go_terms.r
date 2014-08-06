# Import libraries
library(lsr)
library(matrixStats)

# Set working directory
wd  <- '~/Documents/projects_data_analysis/dream/' 

setwd(wd)
source('emanuel/dream2014_Import_Datasets.R')

samples <- colnames(exp.train.raw)
genes <- intersect(rownames(exp.train.raw), rownames(achilles.training))

# Read GO terms
go.terms <- read.table(file='emanuel/go_terms/c2.all.v4.0.symbols.gmt.tab', sep='\t', stringsAsFactors=F, check.names=F, quote='', header=T)
go.terms <- go.terms[which(go.terms$gene %in% genes),]

sample <- samples[3]

# Normalise sample gex
sample.gex <- data.frame(gene=genes, exp=exp.train.raw[genes, sample])
sample.gex <- sample.gex[order(sample.gex$exp, decreasing=F),]

sample.gex.norm.fun <- ecdf(sample.gex$exp)
sample.gex.norm <- data.frame(row.names=sample.gex$gene, exp=sample.gex.norm.fun(sample.gex$exp))

# Normalise sample essentiality
sample.ess <- data.frame(gene=genes, exp=achilles.training[genes ,sample])
sample.ess <- sample.ess[order(sample.ess$exp, decreasing=T),]

sample.ess.norm.fun <- ecdf(sample.ess$exp)
sample.ess.norm <- data.frame(row.names=sample.ess$gene, exp=sample.ess.norm.fun(sample.ess$exp))

term.gsea <- do.call(rbind, lapply(unique(go.terms$term), function(term) {
  term.genes <- go.terms[which(go.terms$term == term),]$gene
  es.ess <- wGSEA(rownames(sample.gex.norm), sample.gex.norm$exp, term.genes)
  es.exp <- wGSEA(rownames(sample.ess.norm), sample.ess.norm$exp, term.genes)
  data.frame(row.names=term, ess=es.ess, exp=es.exp)
}))