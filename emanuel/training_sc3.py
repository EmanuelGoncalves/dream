__author__ = 'emanuel'

import numpy as np
import operator
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV, PassiveAggressiveRegressor
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.cross_validation import ShuffleSplit
from sklearn.metrics import make_scorer
from sklearn.svm import *
from sklearn.neural_network import BernoulliRBM
from pandas import DataFrame, Series
from dream_2014_functions import read_data_sets, leader_board_cell_lines, training_cell_lines


def spearm_cor_func(expected, pred):
    return spearmanr(expected, pred)[0]

# Import data
exp, cnv, ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Split training data-set in two
train_exp = exp.loc[training_cell_lines, ]
train_cnv = cnv.loc[training_cell_lines, ]
train_ess = ess.loc[training_cell_lines, ]

pred_exp = exp.loc[leader_board_cell_lines, ]
pred_cnv = cnv.loc[leader_board_cell_lines, ]
pred_ess = ess.loc[leader_board_cell_lines, ].T

# Predicted genes
samples = pred_ess.axes[1]

# Configuration
X_train_pre = train_exp
X_test_pre = pred_exp

# Filter by coeficient variation
var_thres = VarianceThreshold(.625).fit(X_train_pre)
X_train_pre = X_train_pre.loc[:, var_thres.get_support()]
X_test_pre = X_test_pre.loc[:, var_thres.get_support()]

# Filter by correlation
features = X_train_pre.columns.values

train_exp_corcoef = np.corrcoef(X_train_pre.T)
train_exp_corcoef = np.triu(train_exp_corcoef)
train_exp_corcoef = np.where(train_exp_corcoef > 0.80)

train_exp_corcoef = [(train_exp_corcoef[0][i], train_exp_corcoef[1][i]) for i in range(len(train_exp_corcoef[0])) if train_exp_corcoef[0][i] != train_exp_corcoef[1][i]]
features_to_remove = set(X_train_pre.columns[x[1]] for x in train_exp_corcoef)

X_train_pre = X_train_pre.drop(features_to_remove, 1)
X_test_pre = X_test_pre.drop(features_to_remove, 1)

# Prepare features
features = X_train_pre.columns
important_features = []

for gene in prioritized_genes:
    # Assemble prediction variables
    X_train = X_train_pre
    y_train = train_ess.ix[:, gene]
    X_test = X_test_pre

    # Feature selection
    fs = SelectKBest(f_regression, k=150)
    X_train = fs.fit_transform(X_train, y_train)
    X_test = fs.transform(X_test)
    gene_features = features[fs.get_support()]

    # Store gene features
    important_features.append(gene_features.values)

    print gene, X_test.shape

# Filter features
important_features = Series([feature for feature_list in important_features for feature in feature_list])
important_features_top_100 = [x for x in important_features.value_counts().head(100).index]

predictions = DataFrame(None, index=prioritized_genes, columns=samples)
spearman = make_scorer(spearm_cor_func, greater_is_better=True)

# Assemble prediction variables
X_train = X_train_pre.loc[:, important_features_top_100]
X_test = X_test_pre.loc[:, important_features_top_100]

for gene in prioritized_genes:
    y_train = train_ess.ix[:, gene]

    y_pred = np.mean([
        PassiveAggressiveRegressor(epsilon=0.01, n_iter=3).fit(X_train, y_train).predict(X_test),
        RidgeCV(normalize=True).fit(X_train, y_train).predict(X_test),
    ], axis=0)

    print gene, X_train.shape

    # Store results
    predictions.ix[gene] = y_pred


# Calculate score
correlations = []
for gene in prioritized_genes:
    correlations.append(spearmanr(predictions.loc[gene], pred_ess.loc[gene])[0])

# Register run result
score = '%.5f' % np.mean(correlations)

print score