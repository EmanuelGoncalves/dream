__author__ = 'emanuel'

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.cluster.vq import kmeans2
from scipy.stats import boxcox
from statsmodels.distributions import ECDF as ecdf
from sklearn.linear_model import *
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA, NMF, RandomizedPCA
from sklearn.random_projection import johnson_lindenstrauss_min_dim, GaussianRandomProjection
from sklearn.lda import LDA
from sklearn.cluster import FeatureAgglomeration
from sklearn.grid_search import GridSearchCV
from sklearn.pipeline import Pipeline, FeatureUnion
from sklearn.feature_selection import *
from sklearn.isotonic import IsotonicRegression
from sklearn.ensemble import *
from sklearn.tree import *
from sklearn.neural_network import *
from sklearn.tree import DecisionTreeRegressor
from sklearn.cross_decomposition import PLSRegression
from sklearn.manifold import *
from sklearn.svm import *
from sklearn.cross_decomposition import PLSCanonical, PLSRegression, CCA
from sklearn.cross_validation import *
from sklearn.metrics import *
from sklearn.neighbors import *
from sklearn.feature_selection import f_regression, SelectKBest, chi2
from pandas import DataFrame, Series
from dream_2014_functions import read_data_sets, leader_board_cell_lines, training_cell_lines


def hill_function(matrix, hill_coef=7):
    return 1 / ((np.median(matrix, axis=0) / matrix) ** hill_coef + 1)


def spearm_cor_func(y_true, y_pred):
    return spearmanr(y_true, y_pred)[0]


def discretise_cnv(matrix, filter_sd=True, lower_bound=-1, upper_bound=1):
    matrix_discrete = DataFrame(0, index=matrix.axes[0], columns=matrix.axes[1])
    matrix_discrete[matrix <= lower_bound] = -1.2
    matrix_discrete[matrix >= upper_bound] = 1.2
    return matrix_discrete.loc[:, matrix_discrete.std() != 0] if filter_sd else matrix_discrete


def cor_exp_ess(exp, ess):
    cor = DataFrame(np.nan, index=ess.columns, columns=['cor', 'pvalue'])

    for gene in ess.columns:
        if gene in exp.columns:
            cor.loc[gene] = spearmanr(ess[gene], exp[gene])

    return cor.dropna()


def ecdf_norm(matrix):
    norm_matrix = DataFrame(0, index=matrix.axes[0], columns=matrix.axes[1])
    for column in matrix.columns:
        norm_matrix.loc[:, column] = ecdf(matrix.loc[:, column])(matrix.loc[:, column])
    return norm_matrix


def discretise_kmeans(matrix, n_cluster=3):
    matrix_discrete = DataFrame(
        [kmeans2(matrix.ix[:, i].values, n_cluster)[1] for i in range(matrix.shape[1])],
        index=matrix.columns,
        columns=matrix.index
    )

    return matrix_discrete.T


def f_regression_aux (X_test, y_train, center=False):
    return f_regression(X_test, y_train, center)

# Import data-sets
exp, cnv, ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Remove columns names with NaNs
exp = exp.loc[:, [str(x) != 'nan' for x in exp.columns]]
leader_exp = leader_exp.loc[:, [str(x) != 'nan' for x in leader_exp.columns]]

# Split training data-set in two
train_exp = exp.loc[training_cell_lines, ]
train_cnv = cnv.loc[training_cell_lines, ]
train_ess = ess.loc[training_cell_lines, ]

pred_exp = exp.loc[leader_board_cell_lines, ]
pred_cnv = cnv.loc[leader_board_cell_lines, ]
pred_ess = ess.loc[leader_board_cell_lines, ].T

# Predicted genes
genes = pred_ess.axes[0]
# genes = prioritized_genes
samples = pred_ess.axes[1]

# Configurations
predictions = DataFrame(None, index=genes, columns=samples)
spearman = make_scorer(spearm_cor_func, greater_is_better=True)

X_train_pre = train_exp
X_test_pre = pred_exp

best_epsilon = 0.01
best_n_iter = 3
best_k = 3500

# Filter by coeficient variation
# features_to_keep = X_train_pre.std() / X_train_pre.mean() > 0.1
# X_train_pre = X_train_pre.loc[:, features_to_keep.values]
# X_test_pre = X_test_pre.loc[:, features_to_keep.values]

var_thres = VarianceThreshold(.625).fit(train_exp)
train_exp = train_exp.loc[:, var_thres.get_support()]
pred_exp = pred_exp.loc[:, var_thres.get_support()]

var_thres = VarianceThreshold(.2).fit(train_cnv)
train_cnv = train_cnv.loc[:, var_thres.get_support()]
pred_cnv = pred_cnv.loc[:, var_thres.get_support()]

# X_train_pre = train_exp.join(train_cnv, rsuffix='_cnv')
# X_test_pre = pred_exp.join(pred_cnv, rsuffix='_cnv')

X_train_pre = train_exp
X_test_pre = pred_exp

# Filter by correlation
train_exp_corcoef = np.corrcoef(X_train_pre.T)
train_exp_corcoef = np.triu(train_exp_corcoef)
train_exp_corcoef = np.where(train_exp_corcoef > 0.85)

train_exp_corcoef = [(train_exp_corcoef[0][i], train_exp_corcoef[1][i]) for i in range(len(train_exp_corcoef[0])) if train_exp_corcoef[0][i] != train_exp_corcoef[1][i]]
features_to_remove = set(X_train_pre.columns[x[1]] for x in train_exp_corcoef)

X_train_pre = X_train_pre.drop(features_to_remove, 1)
X_test_pre = X_test_pre.drop(features_to_remove, 1)

corrs = []

for gene in genes:
    # Assemble prediction variables
    X_train = X_train_pre
    y_train = train_ess.ix[:, gene].values
    X_test = X_test_pre

    # Feature selection
    fs = SelectKBest(f_regression, k=best_k).fit(X_train, y_train)
    X_train = X_train.loc[:, fs.get_support()]
    X_test = X_test.loc[:, fs.get_support()]

    y_pred = np.mean([
        PassiveAggressiveRegressor(epsilon=0.01, n_iter=3).fit(X_train, y_train).predict(X_test),
        RidgeCV(normalize=True).fit(X_train, y_train).predict(X_test),
        NuSVR().fit(X_train, y_train).predict(X_test),
    ], axis=0)

    # y_pred = PassiveAggressiveRegressor(epsilon=best_epsilon, n_iter=best_n_iter).fit(X_train, y_train).predict(X_test)
    # y_pred = RidgeCV(normalize=True).fit(X_train, y_train).predict(X_test)

    predictions.loc[gene, :] = y_pred

    cor = spearmanr(predictions.loc[gene, :], pred_ess.loc[gene])[0]
    corrs.append(cor)
    print gene, X_train.shape, cor, '\t\t', np.mean(corrs)

# Calculate score
correlations = [spearmanr(predictions.loc[gene], pred_ess.loc[gene])[0] for gene in genes]

# Register run result
score = '%.5f' % np.mean(correlations)
print score