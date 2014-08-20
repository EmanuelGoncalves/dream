library(pheatmap)

wd  <- '~/Documents/projects_data_analysis/dream/' 

# Set working directory
source(paste(wd, 'dream2014-functions.R', sep=''))

# Read training gene expression
exp.train.raw <- read.gct.file(paste(wd, '_data_phase2/CCLE_expression_training_phase2.gct', sep=''))
exp.train.raw <- exp.train.raw[which(exp.train.raw$Description != 'NaN'),]
probes.to.remove <- c('55016_at', '64757_at', '51257_at', '54996_at')
exp.train.raw <- exp.train.raw[-which(exp.train.raw$Name %in% probes.to.remove),]
exp.train.raw <- exp.train.raw[,-1]
row.names(exp.train.raw) <- exp.train.raw$Description
exp.train.raw <- exp.train.raw[,-1]
colnames(exp.train.raw)[1] <- '143B'
colnames(exp.train.raw)[2] <- '769P'

# Read leader gene expression
exp.leader.raw <- read.gct.file(paste(wd, '_data_phase2/CCLE_expression_leaderboard_phase2.gct', sep=''))
exp.leader.raw <- exp.leader.raw[which(exp.leader.raw$Description != 'NaN'),]
probes.to.remove <- c('55016_at', '64757_at', '51257_at', '54996_at')
exp.leader.raw <- exp.leader.raw[-which(exp.leader.raw$Name %in% probes.to.remove),]
exp.leader.raw <- exp.leader.raw[,-1]
row.names(exp.leader.raw) <- exp.leader.raw$Description
exp.leader.raw <- exp.leader.raw[,-1]

# Read annotation files
training.annot.file <- read.table(file=paste(wd, '_data_phase2/Cell_line_annotation_training_phase2.txt', sep=''), sep='\t', header=T, row.names=1)
training.annot.file$dream <- 'training'

leader.annot.file <- read.table(file=paste(wd, '_data_phase2/Cell_line_annotation_leaderboard_phase2.txt', sep=''), sep='\t', header=T, row.names=1)
leader.annot.file$dream <- 'leader'

annot <- rbind(training.annot.file, leader.annot.file)

# Assemble gene expression matrix
exp <- cbind(exp.train.raw, exp.leader.raw)

exp.norm <- t(apply(exp, 1, scale))
colnames(exp.norm) <- colnames(exp)

annot.ph <- annot[,c(-1)]

pheatmap(exp.norm, show_rownames=F, cluster_rows=F, annotation=annot.ph)

exp.norm.d <- dist(t(as.matrix(exp.norm)), method = 'euclidean')
exp.norm.fit <- hclust(exp.norm.d)

plot(exp.norm.fit)
rect.hclust(exp.norm.fit, k = 6, border = 'blue')
