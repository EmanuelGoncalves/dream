source("dream2014-functions.R")

#### READING DATA ####

# checking cell line annotations for the different phases
cell.line.annot.tr.1 = read.annotations("data/phase1/Cell_line_annotation_training.txt")
dim(cell.line.annot.tr.1)
names(cell.line.annot.tr.1)
names.lines.tr.1 = cell.line.annot.tr.1$Name

cell.line.annot.ts.1 = read.annotations("data/phase1/Cell_line_annotation_leaderboard.txt")
dim(cell.line.annot.ts.1)
names(cell.line.annot.ts.1)
names.lines.ts.1 = cell.line.annot.ts.1$Name

cell.line.annot.tr.2 = read.annotations("data/phase2/Cell_line_annotation_training_phase2.txt")
dim(cell.line.annot.tr.2)
names.lines.tr.2 = cell.line.annot.tr.2$Name
names.lines.tr.2[5] = toupper(names.lines.tr.2[5]) # seems to be an error

shared.tr = intersect(names.lines.tr.1, names.lines.tr.2)
length(shared.tr)
# 45

cell.line.annot.ts.2 = read.annotations("data/phase2/Cell_line_annotation_leaderboard_phase2.txt")
dim(cell.line.annot.ts.2)
names.lines.ts.2 = cell.line.annot.ts.2$Name

shared.ts = intersect(names.lines.ts.1, names.lines.ts.2)
length(shared.ts)
# 22

cell.line.annot.tr.3 = read.annotations("data/phase3/Cell_line_annotation_training_phase3.txt")
dim(cell.line.annot.tr.3)
names.lines.tr.3 = cell.line.annot.tr.3$Name

shared.tr.3 = intersect(names.lines.tr.3, names.lines.tr.2)
length(shared.tr.3)
# 66

length(intersect(new.tr.3, names.lines.ts.2))
# 33

cell.line.annot.ts.3 = read.annotations("data/phase3/Cell_line_annotation_finaltest_phase3.txt")
dim(cell.line.annot.ts.3)
names.lines.ts.3 = cell.line.annot.ts.3$Name



# Loading expression data

exp.train.raw = read.gct.file("data/CCLE_expression_training.gct")
dim(exp.train.raw)
#18960  47

exp.leader.raw = read.gct.file("data/CCLE_expression_leaderboard.gct")
dim(exp.leader.raw)
# 18960 24

sum(duplicated(exp.train.raw$Description))
# 299 genes have more than one probe

exp.annot = read.annotations("data/CCLE_expression_gene_annotations.txt")
dim(exp.annot)
# 18960     2

which(exp.train.raw$Description != exp.annot$Gene_Name)
# 26 are different ... keeping the ones on Description for now ...

range(exp.train.raw[,3:47])
# 2.7766 15.1182
boxplot(exp.train.raw[,3:47])

# Loading copy number data

cn.train.raw = read.gct.file("data/CCLE_copynumber_training.gct")
dim(cn.train.raw)
# [1] 23288    47
head(cn.train.raw)
which(cn.train.raw$Name != cn.train.raw$Description)
# none !
sum(duplicated(cn.train.raw$Name))
# no duplicates

cn.leader.raw = read.gct.file("data/CCLE_copynumber_leaderboard.gct")
dim(cn.leader.raw)
# 23288 24

cn.annot = read.annotations("data/CCLE_copynumber_gene_annotations.txt")
dim(cn.annot)
# 21217     4
head(cn.annot)

sum(duplicated(cn.annot$Gene_Name))
# 0

range(cn.train.raw[,3:47])
# -8.2165  4.6388
boxplot(cn.train.raw[,3:47])

# Loading scores

scores.train.raw = read.gct.file("data/Achilles_v2.9_training.gct")
dim(scores.train.raw)
#[1] 14557    47
head(scores.train.raw)
which(scores.train.raw$Name != scores.train.raw$Description)
# all equal

range(scores.train.raw[,3:47])
boxplot(scores.train.raw[,3:47])

# list of genes with score values
genelist = scores.train.raw$Name
length(genelist)
sum(duplicated(genelist))

# list of genes from the ones with score that have expression data
length(which(exp.annot$Gene_Name %in% genelist))
#12974
genesWithoutExp = genelist [which(!(genelist %in% exp.annot$Gene_Name)) ]
genesWithExp = genelist [which(genelist %in% exp.annot$Gene_Name)]
length(genesWithExp)
# 12974
length(genesWithoutExp)
#1583

# list of genes from the ones with score that have copy number data
length(which(cn.annot$Gene_Name %in% genelist))
# 13819
genesWithCN = genelist [which(genelist %in% cn.annot$Gene_Name)]
genesWithoutCN = genelist [which(!(genelist %in% cn.annot$Gene_Name))]
# 738
genesWithoutInfo = intersect(genesWithoutCN, genesWithoutExp)
length(genesWithoutInfo)
# 360

# kNN example results

knn.results.raw = read.gct.file("submissions/knn_k_3_gene_expression_leaderboard_prediction.gct")
dim(knn.results.raw)
# 14557    24

## ML input matrix

ge.tr.mat = t(exp.train.raw[,3:47])
colnames(ge.tr.mat) = exp.train.raw[,2]
rownames(ge.tr.mat) = colnames(exp.train.raw)[3:47]

cn.tr.mat = t(cn.train.raw[,3:47])
colnames(cn.tr.mat) = cn.train.raw[,1]
rownames(cn.tr.mat) = colnames(cn.train.raw)[3:47]

input.tr.mat = cbind(ge.tr.mat, cn.tr.mat)
rownames(input.tr.mat) = rownames(ge.tr.mat)
colnames(input.tr.mat) = c(paste("E",colnames(ge.tr.mat), sep=""), paste("CN", colnames(cn.tr.mat), sep=""))
dim(input.tr.mat)

output.tr.mat = t(scores.train.raw[,3:47])
rownames(output.tr.mat) = rownames(ge.tr.mat)
colnames(output.tr.mat) = scores.train.raw[,1]

ge.ts.mat = t(exp.leader.raw[,3:24])
colnames(ge.ts.mat) = exp.leader.raw[,2]
rownames(ge.ts.mat) = colnames(exp.leader.raw)[3:24]

cn.ts.mat = t(cn.leader.raw[,3:24])
colnames(cn.ts.mat) = cn.leader.raw[,1]
rownames(cn.ts.mat) = colnames(cn.leader.raw)[3:24]

input.ts.mat = cbind(ge.ts.mat, cn.ts.mat)
rownames(input.ts.mat) = rownames(ge.ts.mat)
colnames(input.ts.mat) = c(paste("E",colnames(ge.ts.mat), sep=""), paste("CN", colnames(cn.ts.mat), sep=""))
dim(input.ts.mat)

# normalization

norm.input.tr.mat = scale(input.tr.mat)
colnames(norm.input.tr.mat) = NULL
norm.input.ts.mat = scale(input.ts.mat)
colnames(norm.input.ts.mat) = NULL
norm.output.tr.mat = scale(output.tr.mat)


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

 
