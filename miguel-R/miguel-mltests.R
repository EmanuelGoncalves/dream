# normalization

norm.input.tr.mat = scale(input.tr.mat)
colnames(norm.input.tr.mat) = NULL
norm.input.ts.mat = scale(input.ts.mat)
colnames(norm.input.ts.mat) = NULL
norm.output.tr.mat = scale(output.tr.mat)

# kNN example results

knn.results.raw = read.gct.file("submissions/knn_k_3_gene_expression_leaderboard_prediction.gct")
dim(knn.results.raw)
# 14557    24

# experiments with kNN reg (from FNN package)
require(FNN)
res.knn.1 = matrix(NA, nrow(input.ts.mat), ncol(output.tr.mat))
for (g in 1:ncol(output.tr.mat) ) {
  res.knn.gene = knn.reg(input.tr.mat, input.ts.mat, output.tr.mat[,g], k = 1)
  res.knn.1[,g] = res.knn.gene$pred
}

write.res.gct(res.knn.1, colnames(output.tr.mat), rownames(input.ts.mat), outfile = "submissions/res.knn.1.gct")

res.knn.3 = matrix(NA, nrow(input.ts.mat), ncol(output.tr.mat))
for (g in 1:ncol(output.tr.mat) ) {
  res.knn.gene = knn.reg(input.tr.mat, input.ts.mat, output.tr.mat[,g], k = 3)
  res.knn.3[,g] = res.knn.gene$pred
}

write.res.gct(res.knn.3, colnames(output.tr.mat), rownames(input.ts.mat), outfile = "submissions/umebi-subm1.gct")

# own kNN approach - based on correlations (much faster)

correl.mat = cor(t(scale(input.ts.mat)), t(scale(input.tr.mat)), method = "spearman")

k = 3
res.our.knn = matrix(NA, nrow(input.ts.mat), ncol(output.tr.mat))
for(cl in 1:nrow(correl.mat)) {
  bestK = order(correl.mat[cl,], decreasing = T)[1:k]
  res.our.knn[cl,] = colMeans(output.tr.mat[bestK,])
}

write.res.gct(res.our.knn, colnames(output.tr.mat), rownames(input.ts.mat), outfile = "submissions/umebi-subm2.gct")

k = 11
res.our.knn.11 = matrix(NA, nrow(input.ts.mat), ncol(output.tr.mat))
for(cl in 1:nrow(correl.mat)) {
  bestK = order(correl.mat[cl,], decreasing = T)[1:k]
  res.our.knn.11[cl,] = colMeans(output.tr.mat[bestK,])
}

write.res.gct(res.our.knn.11, colnames(output.tr.mat), rownames(input.ts.mat), outfile = "submissions/umebi-subm3.gct")

correl.mat.norm = cor(t(scale(norm.input.ts.mat)), t(scale(norm.input.tr.mat)), method = "spearman")

k = 11
res.our.knn.11.v2 = matrix(NA, nrow(input.ts.mat), ncol(output.tr.mat))
for(cl in 1:nrow(correl.mat.norm)) {
  bestK = order(correl.mat.norm[cl,], decreasing = T)[1:k]
  res.our.knn.11.v2[cl,] = colMeans(output.tr.mat[bestK,])
}

write.res.gct(res.our.knn.11.v2, colnames(output.tr.mat), rownames(input.ts.mat), outfile = "submissions/umebi-subm3-2.gct")


# trying SVMs 

require(e1071)
res.svm = matrix(NA, nrow(input.ts.mat), ncol(output.tr.mat))
nf = 1000
for (g in 1:ncol(output.tr.mat) ) {
  correl.feat.scores = cor(norm.input.tr.mat, norm.output.tr.mat[,g])
  bestFeatures = order(correl.feat.scores, decreasing = T)[1:nf]
  tr.set = data.frame(norm.input.tr.mat[,bestFeatures], out = output.tr.mat[,g])
  model.svm = svm(out ~., data = tr.set)
  res.svm[,g] = predict(model.svm, norm.input.ts.mat[,bestFeatures])
  print (g)
}

write.res.gct(res.svm, colnames(output.tr.mat), rownames(input.ts.mat), outfile = "submissions/umebi-subm4-2.gct")

res.svm2 = matrix(NA, nrow(input.ts.mat), ncol(output.tr.mat))
nf = 500
for (g in 1:ncol(output.tr.mat) ) {
  correl.feat.scores2 = abs(cor(norm.input.tr.mat, norm.output.tr.mat[,g]))
  bestFeatures2 = order(correl.feat.scores2, decreasing = T)[1:nf]
  tr.set = data.frame(norm.input.tr.mat[,bestFeatures2], out = output.tr.mat[,g])
  model.svm = svm(out ~., data = tr.set)
  res.svm[,g] = predict(model.svm, norm.input.ts.mat[,bestFeatures2])
  print (g)
}

write.res.gct(res.svm, colnames(output.tr.mat), rownames(input.ts.mat), outfile = "submissions/umebi-subm4-3.gct")

# trying PLS

require(pls)
res.pls = matrix(NA, nrow(input.ts.mat), ncol(output.tr.mat))
nf = 500
for (g in 1:ncol(output.tr.mat) ) {
  correl.feat.scores = cor(norm.input.tr.mat, norm.output.tr.mat[,g])
  bestFeatures = order(correl.feat.scores, decreasing = T)[1:nf]
  tr.set = data.frame(norm.input.tr.mat[,bestFeatures], out = output.tr.mat[,g])
  model.pls = plsr(out ~., data = tr.set, ncomp = 20)
  res.pls[,g] = predict(model.pls, norm.input.ts.mat[,bestFeatures])[,1,20]
  print (g)
}

write.res.gct(res.pls, colnames(output.tr.mat), rownames(input.ts.mat), outfile = "submissions/umebi-subm5.gct")

# trying lasso

require(glmnet)
res.lasso = matrix(NA, nrow(input.ts.mat), ncol(output.tr.mat))
nf = 500
for (g in 1:ncol(output.tr.mat) ) {
  correl.feat.scores = abs(cor(norm.input.tr.mat, norm.output.tr.mat[,g]))
  bestFeatures = order(correl.feat.scores, decreasing = T)[1:nf]
  #tr.set = data.frame(norm.input.tr.mat[,bestFeatures], out = output.tr.mat[,g])
  model.lasso = glmnet(norm.input.tr.mat[,bestFeatures], output.tr.mat[,g], lambda = c(0.01))
  res.lasso[,g] = predict(model.lasso, norm.input.ts.mat[,bestFeatures])
  print (g)
}

write.res.gct(res.lasso, colnames(output.tr.mat), rownames(input.ts.mat), outfile = "submissions/umebi-subm6.gct")

