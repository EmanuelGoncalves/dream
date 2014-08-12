__author__ = 'emanuel'

# Set-up workspace
import sys
sys.path.append('/Users/emanuel/Documents/projects_data_analysis/dream/emanuel/')

import numpy as np
import operator
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from pandas import DataFrame
from dream_2014_functions import read_data_sets, read_features

# Import data
exp, cnv, ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Split training data-set in two
train_exp = exp.iloc[range(0, 50), ]
train_cnv = cnv.iloc[range(0, 50), ]
train_ess = ess.iloc[range(0, 50), ]

pred_exp = exp.iloc[range(51, 66), ]
pred_cnv = cnv.iloc[range(51, 66), ]
pred_ess = ess.iloc[range(51, 66), ].T

# Predicted genes
samples = pred_ess.axes[1]

# Configurations
predictions = DataFrame(None, index=prioritized_genes, columns=samples)
var_fs_thres = 0.25
cv_n = 5

X_train_pre = train_exp
X_test_pre = pred_exp

# Feature selection
features_dict = read_features('submissions/sc2_emanuel_3.txt')
features_soreted = sorted(features_dict.iteritems(), key=operator.itemgetter(1), reverse=True)
features = [features_soreted[i][0] for i in range(100)]

for gene in prioritized_genes:
    # Assemble prediction variables
    X_train = X_train_pre[features]
    y_train = train_ess.ix[:, gene]
    X_test = X_test_pre[features]

    # Normalization
    scaler = StandardScaler().fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)

    # Estimation
    clf = RidgeCV(gcv_mode='auto')
    y_test_pred = clf.fit(X_train, y_train).predict(X_test)

    print gene, X_train.shape, clf.alpha_

    # Store results
    predictions.ix[gene] = y_test_pred

# Calculate score
correlations = []
for gene in prioritized_genes:
    correlations.append(spearmanr(predictions.loc[gene], pred_ess.loc[gene])[0])

# Register run result
score = '%.5f' % np.mean(correlations)

print score