# Import libraries
library(lsr)

# Set working directory
setwd('~/Documents/projects_data_analysis/dream/')
source("dream2014-functions.R")
source("emanuel/wGSEA.R")

# Import data-sets
exp.train.raw <- read.gct.file('data/CCLE_expression_training.gct')
cn.train.raw <- read.gct.file('data/CCLE_copynumber_training.gct')

achilles.training <- read.gct.file('data/Achilles_v2.9_training.gct')

gene.ontology <- read.table(file='emanuel/gene_association.goa_ref_human.txt', sep='\t', stringsAsFactors=F, check.names=F, quote='')

biological.processes <- read.table(file='emanuel/biological_processes_go_terms.txt', sep='\t', stringsAsFactors=F, check.names=F, quote='')
biological.processes <- unique(biological.processes$V1)

molecular.functions <- read.table(file='emanuel/molecular_function_go_terms.txt', sep='\t', stringsAsFactors=F, check.names=F, quote='')
molecular.functions <- unique(molecular.functions$V1)

biological.processes <- c(molecular.functions, biological.processes)

cell.lines.annot <- read.table(file='emanuel/Cell_line_annotation_training.txt', sep='\t', stringsAsFactors=F, check.names=F, quote='', header=T)

# Preprocess exp data-set
exp.train.raw <- exp.train.raw[which(exp.train.raw$Description != 'NaN'),]
probes.to.remove <- c('55016_at', '64757_at', '51257_at', '54996_at')
exp.train.raw <- exp.train.raw[-which(exp.train.raw$Name %in% probes.to.remove),]
exp.train.raw <- exp.train.raw[,-1]
row.names(exp.train.raw) <- exp.train.raw$Description
exp.train.raw <- exp.train.raw[,-1]
colnames(exp.train.raw)[1] = '143B'

# Preprocess gene esssentiality
achilles.training <- achilles.training[,-1]
rownames(achilles.training) <- achilles.training$Description
achilles.training <- achilles.training[,-1]
colnames(achilles.training)[1] = '143B'

# Define predictions
tissues <- unique(cell.lines.annot$Site_primary)
genes <- intersect(rownames(achilles.training), rownames(exp.train.raw))

# Preprocess gene ontology
gene.ontology <- gene.ontology[which(gene.ontology$V5 %in% biological.processes),]
gene.ontology <- gene.ontology[which(gene.ontology$V3 %in% genes),]

genes <- intersect(genes, gene.ontology$V3)
go.terms <- unique(gene.ontology$V5)

go.terms.essentiality <- do.call(cbind, lapply(tissues, function (tissue) {
  samples <- cell.lines.annot[which(cell.lines.annot$Site_primary == tissue),]$Name
  
  genes.scores <- do.call(rbind, lapply(genes, function (gene) {
    values <- unlist(achilles.training[which(rownames(achilles.training) == gene), samples])
    score <- median(values)
    data.frame(gene=gene, score=score)
  }))
  row.names(genes.scores) <- genes.scores$gene
  
  genes.scores <- genes.scores[order(genes.scores$score),]
  
  terms.scores <- do.call(rbind, lapply(go.terms, function (term) {
    term.genes <- unique(gene.ontology[which(gene.ontology$V5 == term),]$V3)
    if (length(term.genes) <= 1) {
      score <- 0
    } else {
      score <- wGSEA(genes.scores$gene, genes.scores$score, term.genes)
    }
    score
  }))
  
  rownames(terms.scores) <- go.terms
  terms.scores <- terms.scores[which(terms.scores != 0),]
  write.table(terms.scores, file=paste('emanuel/tumour_gsea/', tissue, '.csv', sep=''), sep=',', col.names=T, row.names=T, quote=F)
  terms.scores
}))
colnames(go.terms.essentiality) <- tissues
write.table(go.terms.essentiality, file=paste('emanuel/tumour_gsea/all.csv', sep=''), sep=',', col.names=T, row.names=T, quote=F)

# Predictions trainning
samples <- colnames(exp.train.raw)
gene.predictions <- do.call(cbind, lapply(samples, function (sample) {
  tissue <- cell.lines.annot[which(cell.lines.annot$Name == sample),]$Site_primary
  gene.go.terms <- go.terms.essentiality[,tissue]
  
  message(paste(tissue, sample, sep=' - '))
  #dist.fun <- ecdf(exp.train.raw[, sample])
  do.call(rbind, lapply(genes, function (gene) {
    gene.go.terms <- unique(gene.ontology[which(gene.ontology$V3 == gene),]$V5)
    gene.go.terms <- intersect(gene.go.terms, rownames(go.terms.essentiality))
    gene.go.terms <- go.terms.essentiality[gene.go.terms, tissue]
    #gene.exp <- dist.fun(exp.train.raw[gene, sample])
    
    max(gene.go.terms)
  }))
}))
colnames(gene.predictions) <- samples
rownames(gene.predictions) <- genes

predictions.cor <- unlist(lapply(rownames(gene.predictions), function (gene) {
  cor(unlist(achilles.training[gene, colnames(gene.predictions)]), unlist(gene.predictions[gene, colnames(gene.predictions)])*-1, method = 'spearman')
}))

mean(na.omit(predictions.cor))

# Predictions leaderboard
exp.leader.raw <- read.gct.file("data/CCLE_expression_leaderboard.gct")
cell.lines.annot.lb <- read.table(file='data/Cell_line_annotation_leaderboard.txt', sep='\t', stringsAsFactors=F, check.names=F, quote='', header=T)

tissues.lb <- unique(cell.lines.annot.lb$Site_primary)
samples.lb <- cell.lines.annot.lb$Name
genes <- rownames(achilles.training)

gene.predictions.lb <- do.call(cbind, lapply(samples.lb, function (sample) {
  tissue <- cell.lines.annot.lb[which(cell.lines.annot.lb$Name == sample),]$Site_primary
  
  if (tissue %in% colnames(go.terms.essentiality)) {
    gene.go.terms.values <- go.terms.essentiality[,tissue]  
  } else {
    gene.go.terms.values <- rowMedians(go.terms.essentiality)
    names(gene.go.terms.values) <- rownames(go.terms.essentiality)
  }
  
  message(paste(tissue, sample, sep=' - '))
  
  do.call(rbind, lapply(genes, function (gene) {
    gene.go.terms <- unique(gene.ontology[which(gene.ontology$V3 == gene),]$V5)
    
    if (length(gene.go.terms) == 0) {
      as.numeric(median(gene.go.terms.values))
    } else {
      gene.go.terms <- intersect(gene.go.terms, rownames(go.terms.essentiality))
      gene.go.terms <- gene.go.terms.values[gene.go.terms]
      as.numeric(max(gene.go.terms))
    }
  }))
}))
colnames(gene.predictions.lb) <- samples.lb
rownames(gene.predictions.lb) <- genes

gene.predictions.lb <- gene.predictions.lb * -1
gene.predictions.lb <- signif(gene.predictions.lb, digits=4)

gene.predictions.lb[which(is.na(as.numeric(rowSds(gene.predictions.lb)))),] <- colMedians(gene.predictions.lb)

write.res.gct(gene.predictions.lb, rownames(gene.predictions.lb), colnames(gene.predictions.lb), outfile = "emanuel/submissions/umebi_emanuel_1.gct")

predictions.cor <- unlist(lapply(rownames(gene.predictions.lb), function (gene) {
  cor(unlist(gene.predictions.lb[gene, colnames(gene.predictions.lb)]), unlist(gene.predictions.lb[gene, colnames(gene.predictions.lb)]), method = 'spearman')
}))

mean(na.omit(predictions.cor))