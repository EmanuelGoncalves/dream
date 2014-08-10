__author__ = 'emanuel'

# Set-up workspace
import sys
sys.path.append('/Users/emanuel/Documents/projects_data_analysis/dream/emanuel/')

import numpy as np
import operator
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV, ElasticNetCV, LassoCV
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_regression
from pandas import DataFrame
from dream_2014_functions import read_gct, read_features


def register_trainning(method, features='NA', normalization='NA', feature_selection='NA', feature_selection_thres='NA', score='NA', others='NA', sep='\t'):
    with open('emanuel/training.log', 'a') as f:
        f.write(method + sep + features + sep + normalization + sep + feature_selection + sep + feature_selection_thres + sep + score + sep + others + '\n')

# Data-sets files
train_exp_file = 'data/CCLE_expression_training.gct'
train_cnv_file = 'data/CCLE_copynumber_training.gct'
train_ess_file = 'data/Achilles_v2.9_training.gct'
gene_list = 'data/prioritized_gene_list.txt'

# Import data
exp = read_gct(train_exp_file)
cnv = read_gct(train_cnv_file)
ess = read_gct(train_ess_file)

# Split training data-set in two
train_exp = exp.iloc[range(0, 25), ]
train_cnv = cnv.iloc[range(0, 25), ]
train_ess = ess.iloc[range(0, 25), ]

pred_exp = exp.iloc[range(26, 45), ]
pred_cnv = cnv.iloc[range(26, 45), ]
pred_ess = ess.iloc[range(26, 45), ].T

# Predicted genes
genes = np.genfromtxt(gene_list, dtype='str')
samples = pred_ess.axes[1]

# Configurations
predictions = DataFrame(None, index=genes, columns=samples)
var_fs_thres = 0.25
cv_n = 5

X_train_pre = train_exp
X_test_pre = pred_exp

# Feature selection
features_dict = read_features('submissions/sc2_emanuel_3.txt')
features_soreted = sorted(features_dict.iteritems(), key=operator.itemgetter(1), reverse=True)
features = [features_soreted[i][0] for i in range(100)]

for gene in genes:
    # Assemble prediction variables
    X_train = X_train_pre[features]
    y_train = train_ess.ix[:, gene]
    X_test = X_test_pre[features]

    # Normalization
    scaler = StandardScaler().fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)

    # Estimation
    #clf = RidgeCV(gcv_mode='auto')
    #clf = LassoCV()
    clf = ElasticNetCV(l1_ratio=[.1, .5, .7, .9, .95, .99, 1])
    y_test_pred = clf.fit(X_train, y_train).predict(X_test)

    print gene, X_train.shape, clf.alpha_

    # Store results
    predictions.ix[gene] = y_test_pred

# Calculate score
correlations = []
for gene in genes:
    correlations.append(spearmanr(predictions.loc[gene], pred_ess.loc[gene])[0])

# Register run result
score = '%.5f' % np.mean(correlations)

print score