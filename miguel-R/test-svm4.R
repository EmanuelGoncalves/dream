source("miguel-R/dream-io.R")
source("miguel-R/norm.R")

library(e1071)

sink(file = "res-svm.txt")

nf = 1000
res.svm.4 = matrix(NA, nrow(input.ts.mat.2), ncol(output.tr.mat.2))
for (g in 1:ncol(output.tr.mat.2) ) {
  correl.feat.scores2 = abs(cor(norm.input.tr.mat.2, norm.output.tr.mat.2[,g]))
  bestFeatures2 = order(correl.feat.scores2, decreasing = T)[1:nf]
  tr.set = data.frame(norm.input.tr.mat.2[,bestFeatures2], out = output.tr.mat.2[,g])
  model.svm = svm(out ~., data = tr.set)
  res.svm.4[,g] = predict(model.svm, norm.input.ts.mat.2[,bestFeatures2])
  if (g%%100 == 0) print(g)
}

print "NF = 1000"
print obj.function(output.ts.mat.2, res.svm.4)

nf = 2000
res.svm.5 = matrix(NA, nrow(input.ts.mat.2), ncol(output.tr.mat.2))

for (g in 1:ncol(output.tr.mat.2) ) {
  correl.feat.scores2 = abs(cor(norm.input.tr.mat.2, norm.output.tr.mat.2[,g]))
  bestFeatures2 = order(correl.feat.scores2, decreasing = T)[1:nf]
  tr.set = data.frame(norm.input.tr.mat.2[,bestFeatures2], out = output.tr.mat.2[,g])
  model.svm = svm(out ~., data = tr.set)
  res.svm.5[,g] = predict(model.svm, norm.input.ts.mat.2[,bestFeatures2])
  if (g%%100 == 0) print(g)
}

print "NF = 2000"
print obj.function(output.ts.mat.2, res.svm.5)

sink()