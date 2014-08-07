# Import libraries
library(lsr)
library(matrixStats)

# Set working directory
wd  <- '~/Documents/projects_data_analysis/dream/' 
setwd(wd)
source('emanuel/dream2014_Import_Datasets.R')

# Source EN libs
source('EN_demo/lib/LIB_METRICS.R')
source('EN_demo/lib/LIB_ELASTIC_NET.R')
source('EN_demo/lib/LIB_CROSS_VALIDATION.R')

# Set-up
genes <- colnames(achilles.training)
samples <- rownames(achilles.training)

# Assemble features matrices
features <- cbind(exp.train.raw, cn.train.raw)
leader.features <- cbind(exp.leader.raw, cn.leader.raw)

# Filter features matrices
train.feature.col.vars <- colVars(features)
train.feature.quantile <- quantile(train.feature.col.vars, 0.85)
features <- features[,names(which(train.feature.col.vars > train.feature.quantile))]

leader.feature.col.vars <- colVars(leader.features)
leader.feature.quantile <- quantile(leader.feature.col.vars, 0.85)
leader.features <- leader.features[,names(which(leader.feature.col.vars > leader.feature.quantile))]

# Run elastic net and predic leader data-set
leader.preds <- do.call(rbind, lapply(head(genes), function (gene) {
  message(gene)
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
  
  leader.pred <- rowMeans(leader.pred)
}))
rownames(leader.preds) <- genes

# Write predictions
write.res.gct(leader.preds, rownames(leader.pred), colnames(leader.pred), outfile='emanuel/submissions/umebi_emanuel_3.gct')