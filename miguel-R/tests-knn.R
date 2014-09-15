sink("res-knn.txt")
correl.mat = cor(t(scale(norm.input.ts.mat.2)), t(scale(norm.input.tr.mat.2)), 
                 method = "spearman")

k = 1
res.our.knn.1 = matrix(NA, nrow(input.ts.mat.2), ncol(output.tr.mat.2))
for(cl in 1:nrow(res.our.knn.1)) {
  bestK = order(correl.mat[cl,], decreasing = T)[1:k]
  res.our.knn.1[cl,] = colMeans(output.tr.mat.2[bestK,,drop=F])
}

obj.knn.1 = obj.function(output.ts.mat.2, res.our.knn.1)
obj.knn.1
print(obj.knn.1)

k = 3
res.our.knn.3 = matrix(NA, nrow(input.ts.mat.2), ncol(output.tr.mat.2))
for(cl in 1:nrow(res.our.knn.3)) {
  bestK = order(correl.mat[cl,], decreasing = T)[1:k]
  res.our.knn.3[cl,] = colMeans(output.tr.mat.2[bestK,])
}

obj.knn.3 = obj.function(output.ts.mat.2, res.our.knn.3)
obj.knn.3
print(obj.knn.3)

k = 5
res.our.knn.5 = matrix(NA, nrow(input.ts.mat.2), ncol(output.tr.mat.2))
for(cl in 1:nrow(res.our.knn.5)) {
  bestK = order(correl.mat[cl,], decreasing = T)[1:k]
  res.our.knn.5[cl,] = colMeans(output.tr.mat.2[bestK,])
}

obj.knn.5 = obj.function(output.ts.mat.2, res.our.knn.5)
obj.knn.5
print(obj.knn.5)

k = 7
res.our.knn.7 = matrix(NA, nrow(input.ts.mat.2), ncol(output.tr.mat.2))
for(cl in 1:nrow(res.our.knn.7)) {
  bestK = order(correl.mat[cl,], decreasing = T)[1:k]
  res.our.knn.7[cl,] = colMeans(output.tr.mat.2[bestK,])
}

obj.knn.7 = obj.function(output.ts.mat.2, res.our.knn.7)
obj.knn.7
print(obj.knn.7)

k = 9
res.our.knn.9 = matrix(NA, nrow(input.ts.mat.2), ncol(output.tr.mat.2))
for(cl in 1:nrow(res.our.knn.9)) {
  bestK = order(correl.mat[cl,], decreasing = T)[1:k]
  res.our.knn.9[cl,] = colMeans(output.tr.mat.2[bestK,])
}

obj.knn.9 = obj.function(output.ts.mat.2, res.our.knn.9)
obj.knn.9
print(obj.knn.9)

k = 11
res.our.knn.11 = matrix(NA, nrow(input.ts.mat.2), ncol(output.tr.mat.2))
for(cl in 1:nrow(res.our.knn.11)) {
  bestK = order(correl.mat[cl,], decreasing = T)[1:k]
  res.our.knn.11[cl,] = colMeans(output.tr.mat.2[bestK,])
}

obj.knn.11 = obj.function(output.ts.mat.2, res.our.knn.11)
obj.knn.11
print(obj.knn.11)


sink()