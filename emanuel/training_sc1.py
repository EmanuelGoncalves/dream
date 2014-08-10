__author__ = 'emanuel'

# Set-up workspace
import sys
sys.path.append('/Users/emanuel/Documents/projects_data_analysis/dream/emanuel/')

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import PassiveAggressiveRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn import cross_validation
from pandas import DataFrame
from dream_2014_functions import read_gct


def register_trainning(method, features='NA', normalization='NA', feature_selection='NA', feature_selection_thres='NA', score='NA', others='NA', sep='\t'):
    with open('emanuel/training.log', 'a') as f:
        f.write(method + sep + features + sep + normalization + sep + feature_selection + sep + feature_selection_thres + sep + score + sep + others + '\n')


# Data-sets files
train_exp_file = 'data/CCLE_expression_training.gct'
train_cnv_file = 'data/CCLE_copynumber_training.gct'
train_ess_file = 'data/Achilles_v2.9_training.gct'

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
genes = pred_ess.axes[0]
samples = pred_ess.axes[1]

# Configurations
predictions = DataFrame(None, index=genes, columns=samples)
var_fs_thres = 0.25
par_epsilon = 0.01
par_c = 1
cv_n = 5

X_train_pre = train_exp
X_test_pre = pred_exp

var_fs = VarianceThreshold(var_fs_thres)
X_train_pre = var_fs.fit_transform(X_train_pre)
X_test_pre = var_fs.transform(X_test_pre)

fs_methods = 'Variance'
configuration = 'var_fs_thres=' + str(var_fs_thres) + '; ' + 'par_epsilon=' + str(par_epsilon) + '; ' + 'C=' + str(par_c) + '; ' +  'cv_n=' + str(cv_n)
print configuration

for gene in genes:
    # Assemble prediction variables
    X_train = X_train_pre
    y_train = train_ess.ix[:, gene]
    X_test = X_test_pre

    # Normalization
    scaler = StandardScaler().fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)

    print gene, X_train.shape

    # Estimation
    par = PassiveAggressiveRegressor(C=par_c, epsilon=par_epsilon)
    y_test_pred = par.fit(X_train, y_train).predict(X_test)

    # Store results
    predictions.ix[gene] = y_test_pred

# Calculate score
correlations = []
for gene in genes:
    correlations.append(spearmanr(predictions.loc[gene], pred_ess.loc[gene])[0])

# Register run result
score = '%.5f' % np.mean(correlations)

print score
register_trainning(
    'MultiTaskElasticNetCV',
    features='GEX',
    normalization='u=0;s=1',
    feature_selection=fs_methods,
    score=score,
    others=configuration
)