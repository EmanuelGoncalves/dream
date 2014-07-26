# set working directory so that it is pointing to this file 
# (needed for loading) the source the library files.
# setwd("...")

# ----------------------------------------------------------------------
# load demo data
# ----------------------------------------------------------------------
data(iris)
features <- iris[,1:4]
observation <- as.numeric(iris[,5])
names(observation) <- rownames(features)

# ----------------------------------------------------------------------
# source libraries
# ----------------------------------------------------------------------
source("lib/LIB_METRICS.R")
source("lib/LIB_ELASTIC_NET.R")
source("lib/LIB_CROSS_VALIDATION.R")

# ----------------------------------------------------------------------
# get 10-fold cross-validation splits
# ----------------------------------------------------------------------
nFold = 10
cv <- crossvalidation(as.numeric(rownames(iris)), nFold)

# ----------------------------------------------------------------------
# itterate over 10 folds, build 10 models, and predict with
# independent test set, which was never used for training
# ----------------------------------------------------------------------
models <- list()
pred <- c()
obs <- c()
print("train models...")
for (i in 1:length(cv)) {
  print(paste("  -> model", i))
  
  # train model
  input_train <- features[cv[[i]]$train,]
  input_xTrain <- features[cv[[i]]$xTrain,]
  
  obs_train <- observation[cv[[i]]$train]
  obs_xTrain <- observation[cv[[i]]$xTrain]

  models[[paste("model_",i,sep="")]] <- trainEN(input_train, 
                                                obs_train, 
                                                input_xTrain, 
                                                obs_xTrain)
  
  # predict wiht model
  input_test <- features[cv[[i]]$test,]
  pred <- c(pred, predictEN(models[[i]], input_test))
  
  obs <- c(obs, observation[cv[[i]]$test])
}

# ----------------------------------------------------------------------
# prediction versus observation plot
# ----------------------------------------------------------------------
plot(obs, pred)

# ----------------------------------------------------------------------
# plot elastic net with all lambda for first fold
# ----------------------------------------------------------------------
foldX <- 1

plot(models[[foldX]]$model, 
     xvar="lambda", 
     label=TRUE)
varSeqNumber <- sort(rowSums(models[[foldX]]$model$beta == 0), index.return=TRUE)$ix
legend("topright",
       legend=paste("(", varSeqNumber, ") ", rownames(models[[foldX]]$model$beta)[varSeqNumber], sep=""), 
       lty=0,cex=.6)

# ----------------------------------------------------------------------
# plot cross-training
# ----------------------------------------------------------------------
plot(log(models[[foldX]]$model$lambda), models[[foldX]]$train_cIDX, type="l", 
     main="cross-training", xlab="log(lambda)", ylab="c-index")
points(log(models[[foldX]]$model$lambda), models[[foldX]]$train_cIDX)
lines(log(models[[foldX]]$model$lambda), models[[foldX]]$xTrain_cIDX, col="red")
points(log(models[[foldX]]$model$lambda), models[[foldX]]$xTrain_cIDX, col="red")
abline(v=log(models[[foldX]]$lambda), col="purple", lwd=3)
legend("bottomleft", 
       legend=c("training", "cross-training", "chosen lambda"), 
       lty=1, lwd=c(1,1,3), pch=c(1,1,NA),
       col=c("black", "red", "purple"))

# ----------------------------------------------------------------------
# plot feature weigths plot
# ----------------------------------------------------------------------
barplot(sort(models[[foldX]]$featWeights),las=2, cex.names=0.5, cex.main=1, cex.lab=3)
