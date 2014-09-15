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

