__author__ = 'emanuel'

import numpy as np
import re
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from sklearn.linear_model import *
from sklearn.svm import *
from sklearn.pipeline import *
from sklearn.preprocessing import *
from sklearn.feature_selection import *
from sklearn.grid_search import GridSearchCV
from sklearn.cross_validation import *
from sklearn.ensemble import *
from sklearn.tree import *
from sklearn.neighbors import *
from sklearn.neural_network import BernoulliRBM
from sklearn.gaussian_process import *
from sklearn.cross_decomposition import *
from sklearn.metrics import *
from sklearn.isotonic import *
from sklearn.cluster import *
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from pandas import DataFrame, Series, read_csv
from dream_2014_functions import read_data_sets, leader_board_cell_lines, training_cell_lines, read_gene_neighbours


def spearm_cor_func(expected, pred):
    return spearmanr(expected, pred)[0]


def hill_function(matrix, hill_coef=4):
    return 1 / ((np.median(matrix, axis=0) / matrix) ** hill_coef + 1)


def count_outliers(matrix):
    outliers_counts = []
    for i in range(len(matrix.columns)):
        Q1 = np.percentile(matrix.ix[:, i], 25)
        Q3 = np.percentile(matrix.ix[:, i], 75)
        IQR = Q3 - Q1

        outliers_counts.append(sum(matrix.ix[:, i] < (Q1 - 1.5 * IQR)) + sum(matrix.ix[:, i] > (Q3 + 1.5 * IQR)))

    return outliers_counts

# Import data-sets
exp, cnv, ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Split training data-set in two
train_exp = exp.loc[training_cell_lines, ]
train_cnv = cnv.loc[training_cell_lines, ]
train_ess = ess.loc[training_cell_lines, ]

pred_exp = exp.loc[leader_board_cell_lines, ]
pred_cnv = cnv.loc[leader_board_cell_lines, ]
pred_ess = ess.loc[leader_board_cell_lines, ].T

# Configurations
predictions = DataFrame(None, index=prioritized_genes, columns=pred_ess.axes[1])
spearman = make_scorer(spearm_cor_func, greater_is_better=True)
predictions_features = {}

X_train_pre = train_exp
X_test_pre = pred_exp

# Filter by coeficient variation
var_thres = VarianceThreshold(.65).fit(X_train_pre)
X_train_pre = X_train_pre.loc[:, var_thres.get_support()]
X_test_pre = X_test_pre.loc[:, var_thres.get_support()]

# Filter by correlation
# features = X_train_pre.columns.values
#
# train_exp_corcoef = np.corrcoef(X_train_pre.T)
# train_exp_corcoef = np.triu(train_exp_corcoef)
# train_exp_corcoef = np.where(train_exp_corcoef > 0.80)
#
# train_exp_corcoef = [(train_exp_corcoef[0][i], train_exp_corcoef[1][i]) for i in range(len(train_exp_corcoef[0])) if train_exp_corcoef[0][i] != train_exp_corcoef[1][i]]
# features_to_remove = set(X_train_pre.columns[x[1]] for x in train_exp_corcoef)
#
# X_train_pre = X_train_pre.drop(features_to_remove, 1)
# X_test_pre = X_test_pre.drop(features_to_remove, 1)

corrs = []

for gene in prioritized_genes:
    # Assemble prediction variables
    X_train = X_train_pre
    y_train = train_ess.ix[:, gene].values
    X_test = X_test_pre

    fs = SelectKBest(f_regression).fit(X_train, y_train)
    X_train = X_train.loc[:, fs.get_support()]
    X_test = X_test.loc[:, fs.get_support()]

    y_pred = np.mean([
        PassiveAggressiveRegressor(epsilon=0.01, n_iter=3).fit(X_train, y_train).predict(X_test),
        RidgeCV(normalize=True).fit(X_train, y_train).predict(X_test),
    ], axis=0)

    predictions.loc[gene, :] = y_pred

    cor = spearmanr(predictions.loc[gene, :], pred_ess.loc[gene])[0]
    corrs.append(cor)
    print gene, X_train.shape, cor, '\t\t', np.mean(corrs)

# Calculate score
correlations = [spearmanr(predictions.loc[gene], pred_ess.loc[gene])[0] for gene in prioritized_genes]

# Register run result
score = '%.5f' % np.mean(correlations)

print score