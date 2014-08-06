# Set working directory
source(paste(wd, 'emanuel/functions.r', sep=''))
source(paste(wd, 'emanuel/wGSEA.R', sep=''))

# Read training gene expression
exp.train.raw <- read.gct.file(paste(wd, 'data/CCLE_expression_training.gct', sep=''))
exp.train.raw <- exp.train.raw[which(exp.train.raw$Description != 'NaN'),]
probes.to.remove <- c('55016_at', '64757_at', '51257_at', '54996_at')
exp.train.raw <- exp.train.raw[-which(exp.train.raw$Name %in% probes.to.remove),]
exp.train.raw <- exp.train.raw[,-1]
row.names(exp.train.raw) <- exp.train.raw$Description
exp.train.raw <- exp.train.raw[,-1]
colnames(exp.train.raw)[1] = '143B'
exp.train.raw <- t(exp.train.raw)

# Read training copy number variation
cn.train.raw <- read.gct.file(paste(wd, 'data/CCLE_copynumber_training.gct', sep=''))
cn.train.raw <- cn.train.raw[,-1]
row.names(cn.train.raw) <- cn.train.raw$Description
cn.train.raw <- cn.train.raw[,-1]
colnames(cn.train.raw)[1] = '143B'
cn.train.raw <- t(cn.train.raw)

# Read training gene essentiality data-set
achilles.training <- read.gct.file(paste(wd, 'data/Achilles_v2.9_training.gct', sep=''))
achilles.training <- achilles.training[,-1]
rownames(achilles.training) <- achilles.training$Description
achilles.training <- achilles.training[,-1]
colnames(achilles.training)[1] = '143B'
achilles.training <- t(achilles.training)

# Read leader gene expression
exp.leader.raw <- read.gct.file(paste(wd, 'data/CCLE_expression_leaderboard.gct', sep=''))
exp.leader.raw <- exp.leader.raw[which(exp.leader.raw$Description != 'NaN'),]
probes.to.remove <- c('55016_at', '64757_at', '51257_at', '54996_at')
exp.leader.raw <- exp.leader.raw[-which(exp.leader.raw$Name %in% probes.to.remove),]
exp.leader.raw <- exp.leader.raw[,-1]
row.names(exp.leader.raw) <- exp.leader.raw$Description
exp.leader.raw <- exp.leader.raw[,-1]
exp.leader.raw <- t(exp.leader.raw)

# Read leader copy number variation
cn.leader.raw <- read.gct.file('data/CCLE_copynumber_leaderboard.gct')
cn.leader.raw <- cn.leader.raw[,-1]
row.names(cn.leader.raw) <- cn.leader.raw$Description
cn.leader.raw <- cn.leader.raw[,-1]
cn.leader.raw <- t(cn.leader.raw)

# Cell lines annotations
train.annot <- read.table(file=paste(wd, 'data/Cell_line_annotation_training.txt', sep=''), sep='\t', stringsAsFactors=F, check.names=F, quote='', header=T)
leader.annot <- read.table(file=paste(wd, 'data/Cell_line_annotation_leaderboard.txt', sep=''), sep='\t', stringsAsFactors=F, check.names=F, quote='', header=T)

tissues <- unique(c(unique(train.annot$Site_primary), unique(leader.annot$Site_primary)))
tumours <- unique(c(unique(train.annot$Type), unique(leader.annot$Type)))

# Build train tissue type feature matrix
train.tissue.feature <- do.call(rbind, lapply(train.annot$Name, function (sample) {
  sample.tissue <- train.annot[which(train.annot$Name==sample),]$Site_primary
  as.numeric(sample.tissue == tissues)
}))
colnames(train.tissue.feature) <- tissues
rownames(train.tissue.feature) <- train.annot$Name

# Build train tumour type feature matrix
train.tumour.feature <- do.call(rbind, lapply(train.annot$Name, function (sample) {
  sample.tumour <- train.annot[which(train.annot$Name==sample),]$Type
  as.numeric(sample.tumour == tumours)
}))
colnames(train.tumour.feature) <- tumours
rownames(train.tumour.feature) <- train.annot$Name

# Build leader tissue type feature matrix
leader.tissue.feature <- do.call(rbind, lapply(leader.annot$Name, function (sample) {
  sample.tissue <- leader.annot[which(leader.annot$Name==sample),]$Site_primary
  as.numeric(sample.tissue == tissues)
}))
colnames(leader.tissue.feature) <- tissues
rownames(leader.tissue.feature) <- leader.annot$Name

# Build leader tumour type feature matrix
leader.tumour.feature <- do.call(rbind, lapply(leader.annot$Name, function (sample) {
  sample.tumour <- leader.annot[which(leader.annot$Name==sample),]$Type
  as.numeric(sample.tumour == tumours)
}))
colnames(leader.tumour.feature) <- tumours
rownames(leader.tumour.feature) <- leader.annot$Name