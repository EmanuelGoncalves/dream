
norm.input.tr.mat.2 = scale(input.tr.mat.2)
colnames(norm.input.tr.mat.2) = NULL
norm.input.ts.mat.2 = scale(input.ts.mat.2)
colnames(norm.input.ts.mat.2) = NULL
norm.output.tr.mat.2 = scale(output.tr.mat.2)

correl.mat = cor(t(scale(norm.input.ts.mat.2)), t(scale(norm.input.tr.mat.2)), 
                 method = "spearman")

k = 3
res.our.knn = matrix(NA, nrow(input.ts.mat.2), ncol(output.tr.mat.2))
for(cl in 1:nrow(res.our.knn)) {
  bestK = order(correl.mat[cl,], decreasing = T)[1:k]
  res.our.knn[cl,] = colMeans(output.tr.mat.2[bestK,])
}

obj.function(output.ts.mat.2, res.our.knn)
# 0.1176

correl.mat.abs = abs(correl.mat)

res.our.knn.2 = matrix(NA, nrow(input.ts.mat.2), ncol(output.tr.mat.2))
for(cl in 1:nrow(res.our.knn.2)) {
  bestK = order(correl.mat.abs[cl,], decreasing = T)[1:k]
  res.our.knn.2[cl,] = colMeans(output.tr.mat.2[bestK,])
}
obj.function(output.ts.mat.2, res.our.knn.2)
# 0.087

require(e1071)
res.svm = matrix(NA, nrow(input.ts.mat.2), ncol(output.tr.mat.2))
nf = 100
for (g in 1:ncol(output.tr.mat.2) ) {
  correl.feat.scores = cor(norm.input.tr.mat.2, norm.output.tr.mat.2[,g])
  bestFeatures = order(correl.feat.scores, decreasing = T)[1:nf]
  tr.set = data.frame(norm.input.tr.mat.2[,bestFeatures], out = output.tr.mat.2[,g])
  model.svm = svm(out ~., data = tr.set)
  res.svm[,g] = predict(model.svm, norm.input.ts.mat.2[,bestFeatures])
  print (g)
}

obj.function(output.ts.mat.2, res.svm)
# 0.1398

res.svm.2 = matrix(NA, nrow(input.ts.mat.2), ncol(output.tr.mat.2))
nf = 500
for (g in 1:ncol(output.tr.mat.2) ) {
  correl.feat.scores2 = abs(cor(norm.input.tr.mat.2, norm.output.tr.mat.2[,g]))
  bestFeatures2 = order(correl.feat.scores2, decreasing = T)[1:nf]
  tr.set = data.frame(norm.input.tr.mat.2[,bestFeatures2], out = output.tr.mat.2[,g])
  model.svm = svm(out ~., data = tr.set)
  res.svm.2[,g] = predict(model.svm, norm.input.ts.mat.2[,bestFeatures2])
  print (g)
}

obj.function(output.ts.mat.2, res.svm.2)
# 0.1521

res.svm.3 = matrix(NA, nrow(input.ts.mat.2), ncol(output.tr.mat.2))
nf = 1000
for (g in 1:ncol(output.tr.mat.2) ) {
  correl.feat.scores2 = cor(norm.input.tr.mat.2, norm.output.tr.mat.2[,g])
  bestFeatures2 = order(correl.feat.scores2, decreasing = T)[1:nf]
  tr.set = data.frame(norm.input.tr.mat.2[,bestFeatures2], out = output.tr.mat.2[,g])
  model.svm = svm(out ~., data = tr.set)
  res.svm.3[,g] = predict(model.svm, norm.input.ts.mat.2[,bestFeatures2])
  if (g%%100 == 0) print(g)
}

obj.function(output.ts.mat.2, res.svm.3)
# 0.1471


require(pls)
res.pls = matrix(NA, nrow(input.ts.mat.2), ncol(output.tr.mat.2))
nf = 1000
for (g in 1:ncol(output.tr.mat.2) ) {
  correl.feat.scores = abs(cor(norm.input.tr.mat.2, norm.output.tr.mat.2[,g]))
  bestFeatures = order(correl.feat.scores, decreasing = T)[1:nf]
  tr.set = data.frame(norm.input.tr.mat.2[,bestFeatures], out = output.tr.mat.2[,g])
  model.pls = plsr(out ~., data = tr.set, ncomp = 20)
  res.pls[,g] = predict(model.pls, norm.input.ts.mat.2[,bestFeatures])[,1,20]
  print(g)
}

obj.function(output.ts.mat.2, res.pls)
# 0.1379

require(glmnet)
res.lasso = matrix(NA, nrow(input.ts.mat.2), ncol(output.tr.mat.2))
nf = 500
for (g in 1:ncol(output.tr.mat.2) ) {
  correl.feat.scores = abs(cor(norm.input.tr.mat.2, norm.output.tr.mat.2[,g]))
  bestFeatures = order(correl.feat.scores, decreasing = T)[1:nf]
  #tr.set = data.frame(norm.input.tr.mat[,bestFeatures], out = output.tr.mat[,g])
  model.lasso = glmnet(norm.input.tr.mat.2[,bestFeatures], output.tr.mat.2[,g], lambda = c(0.01))
  res.lasso[,g] = predict(model.lasso, norm.input.ts.mat.2[,bestFeatures])
  print (g)
}

obj.function(output.ts.mat.2, res.lasso)
# 0.1271

