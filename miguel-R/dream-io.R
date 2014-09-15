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

length(intersect(names.lines.tr.3, names.lines.ts.2))
# 33

cell.line.annot.ts.3 = read.annotations("data/phase3/Cell_line_annotation_finaltest_phase3.txt")
dim(cell.line.annot.ts.3)
names.lines.ts.3 = cell.line.annot.ts.3$Name


# Loading expression data

exp.train.raw = read.gct.file("data/phase3/CCLE_expression_training_phase3.gct")
dim(exp.train.raw)
#18960  107
colnames(exp.train.raw)[3:4] = sub("X", "", colnames(exp.train.raw[3:4]))

exp.train.raw.1 = exp.train.raw[,names.lines.tr.1]
dim(exp.train.raw.1)
# 18960 45

exp.test.raw.1 = exp.train.raw[,names.lines.ts.1]
dim(exp.test.raw.1)
# 18960 22

exp.train.raw.2 = exp.train.raw[,names.lines.tr.2]
dim(exp.train.raw.2)
# 18960 66

exp.test.raw.2 = exp.train.raw[,names.lines.ts.2]
dim(exp.test.raw.2)
# 18960 33

exp.test.raw = read.gct.file("data/phase3/CCLE_expression_finaltest_phase3.gct")
dim(exp.test.raw)
# 18960 46

exp.annot = read.annotations("data/phase3/CCLE_expression_gene_annotations_phase3.txt")
dim(exp.annot)
# 18960     2


# Loading copy number data

cn.train.raw = read.gct.file("data/phase3/CCLE_copynumber_training_phase3.gct")
dim(cn.train.raw)
# 23288    107
colnames(cn.train.raw)[3:4] = sub("X", "", colnames(cn.train.raw[3:4]))

cn.train.raw.1 = cn.train.raw[,names.lines.tr.1]
dim(cn.train.raw.1)
# 23288 45

cn.train.raw.2 = cn.train.raw[,names.lines.tr.2]
dim(cn.train.raw.2)
# 23288 66

cn.test.raw = read.gct.file("data/phase3/CCLE_copynumber_finaltest_phase3.gct")
dim(cn.test.raw)
# 23288 46

cn.annot = read.annotations("data/phase3/CCLE_copynumber_gene_annotations_phase3.txt")
dim(cn.annot)
# 21217 4

# Loading scores

scores.train.raw = read.gct.file("data/phase3/Achilles_v2.11_training_phase3.gct")
dim(scores.train.raw)
# 14738 107
colnames(scores.train.raw)[3:4] = sub("X", "", colnames(scores.train.raw[3:4]))

scores.train.raw.1 = scores.train.raw[, names.lines.tr.1]
dim(scores.train.raw.1)
# 14738 47

scores.train.raw.2 = scores.train.raw[, names.lines.tr.2]
dim(scores.train.raw.2)
# 14738 66

scores.test.raw.1 = scores.train.raw[, names.lines.ts.1]
dim(scores.test.raw.1)
# 14738 22

scores.test.raw.2 = scores.train.raw[, names.lines.ts.2]
dim(scores.test.raw.2)
# 14738 33

## ML input matrix

ge.tr.mat = t(exp.train.raw[,3:107])
colnames(ge.tr.mat) = exp.train.raw[,2]
rownames(ge.tr.mat) = colnames(exp.train.raw)[3:107]
dim(ge.tr.mat)
# 105 18960

ge.tr.mat.1 = ge.tr.mat[names.lines.tr.1, ]
dim(ge.tr.mat.1)
# 45 18960

ge.tr.mat.2 = ge.tr.mat[names.lines.tr.2, ]
dim(ge.tr.mat.2)
# 66 18960

cn.tr.mat = t(cn.train.raw[,3:107])
colnames(cn.tr.mat) = cn.train.raw[,1]
rownames(cn.tr.mat) = colnames(cn.train.raw)[3:107]

input.tr.mat = cbind(ge.tr.mat, cn.tr.mat)
rownames(input.tr.mat) = rownames(ge.tr.mat)
colnames(input.tr.mat) = c(paste("E",colnames(ge.tr.mat), sep=""), paste("CN", colnames(cn.tr.mat), sep=""))
dim(input.tr.mat)

input.tr.mat.1 = input.tr.mat[names.lines.tr.1, ]
dim(input.tr.mat.1)
# 45 42248

input.tr.mat.2 = input.tr.mat[names.lines.tr.2, ]
dim(input.tr.mat.2)
# 66 42248

input.ts.mat.1 = input.tr.mat[names.lines.ts.1, ]
dim(input.ts.mat.1)
# 22 42248

input.ts.mat.2 = input.tr.mat[names.lines.ts.2, ]
dim(input.ts.mat.2)
# 33 42248

output.tr.mat = t(scores.train.raw[,3:107])
rownames(output.tr.mat) = rownames(ge.tr.mat)
colnames(output.tr.mat) = scores.train.raw[,1]
dim(output.tr.mat)
# 105 14738

output.tr.mat.1 = output.tr.mat[names.lines.tr.1, ]
dim(output.tr.mat.1)
# 45 14738

output.tr.mat.2 = output.tr.mat[names.lines.tr.2, ]
dim(output.tr.mat.2)
# 66 14738

output.ts.mat.1 = output.tr.mat[names.lines.ts.1, ]
dim(output.ts.mat.1)
# 22 14738

output.ts.mat.2 = output.tr.mat[names.lines.ts.2, ]
dim(output.ts.mat.2)
# 33 14738

ge.ts.mat = t(exp.test.raw[,3:46])
colnames(ge.ts.mat) = exp.test.raw[,2]
rownames(ge.ts.mat) = colnames(exp.test.raw)[3:46]

cn.ts.mat = t(cn.test.raw[,3:46])
colnames(cn.ts.mat) = cn.test.raw[,1]
rownames(cn.ts.mat) = colnames(cn.test.raw)[3:46]

input.ts.mat = cbind(ge.ts.mat, cn.ts.mat)
rownames(input.ts.mat) = rownames(ge.ts.mat)
colnames(input.ts.mat) = c(paste("E",colnames(ge.ts.mat), sep=""), paste("CN", colnames(cn.ts.mat), sep=""))
dim(input.ts.mat)
# 44 42248


# OTHER ANALYSES

sum(duplicated(exp.train.raw$Description))
# 299 genes have more than one probe

which(exp.train.raw$Description != exp.annot$Gene_Name)
# 26 are different ...

range(exp.train.raw[,3:47])
# 2.7766 15.1182
boxplot(exp.train.raw[,3:47])

which(cn.train.raw$Name != cn.train.raw$Description)
# none !
sum(duplicated(cn.train.raw$Name))
# no duplicates

range(cn.train.raw[,3:107])
# -8.2165  4.6388
boxplot(cn.train.raw[,3:107])

range(scores.train.raw[,3:107])
boxplot(scores.train.raw[,3:107])

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

