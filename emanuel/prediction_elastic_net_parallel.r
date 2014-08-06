# Read Arguments
args <- commandArgs(trailingOnly = TRUE)
options(echo = TRUE)

wd <- '/nfs/nobackup2/saezgrp/homes/emanuel/dream_2014/'
gene <- args[1]

# Import libraries
library(lsr)
library(matrixStats)

# Set working directory
source(paste(wd, 'dream2014_Import_Datasets.R', sep=''))

# Source EN libs
source(paste(wd, 'EN_demo/lib/LIB_METRICS.R', sep=''))
source(paste(wd, 'EN_demo/lib/LIB_ELASTIC_NET.R', sep=''))
source(paste(wd, 'EN_demo/lib/LIB_CROSS_VALIDATION.R', sep=''))

# Set-up
genes <- colnames(achilles.training)
samples <- rownames(achilles.training)

# Assemble features matrices
features <- cbind(exp.train.raw, cn.train.raw, train.tissue.feature, train.tumour.feature)
leader.features <- cbind(exp.leader.raw, cn.leader.raw, leader.tissue.feature, leader.tumour.feature)

observation <- achilles.training[, gene]

# 10 fold cross-validation
nFold <- 10
cv <- crossvalidation(rownames(features), nFold)

models <- list()
pred <- c()
obs <- c()
for (i in 1:nFold) {
  # train model
  input_train <- features[c(cv[[i]]$train, cv[[i]]$test),]
  input_xTrain <- features[cv[[i]]$xTrain,]
  
  obs_train <- observation[c(cv[[i]]$train, cv[[i]]$test)]
  obs_xTrain <- observation[cv[[i]]$xTrain]
  
  models[[paste('model_',i,sep='')]] <- trainEN(input_train, obs_train, input_xTrain, obs_xTrain)
}

leader.pred <- do.call(cbind, lapply(1:nFold, function (i) {
  predictEN(models[[i]], leader.features)
}))

leader.pred <- rowMedians(leader.pred)

leader.pred <- t(as.matrix(leader.pred))
rownames(leader.pred) <- gene

write.table(leader.pred, file=paste(wd, 'predictions', '/prediction_', gene, '.tab', sep=''), sep='\t', row.names=T, col.names=T, quote=F)