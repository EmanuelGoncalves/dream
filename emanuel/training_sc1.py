__author__ = 'emanuel'

# Set-up workspace
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
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
from sklearn.tree import DecisionTreeRegressor
from sklearn.cross_decomposition import PLSRegression
from sklearn.manifold import *
from sklearn.svm import *
from sklearn.cross_decomposition import PLSCanonical, PLSRegression, CCA
from sklearn.cross_validation import *
from sklearn.metrics import *
from sklearn.feature_selection import f_regression, SelectKBest, chi2
from pandas import DataFrame, Series
from dream_2014_functions import read_data_sets, read_annotations


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

# Import data-sets
exp, cnv, ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Split training data-set in two
train_exp = exp.iloc[range(26, 66), ]
train_cnv = cnv.iloc[range(26, 66), ]
train_ess = ess.iloc[range(26, 66), ]

pred_exp = exp.iloc[range(25), ]
pred_cnv = cnv.iloc[range(25), ]
pred_ess = ess.iloc[range(25), ].T

# Predicted genes
genes = pred_ess.axes[0]
# genes = prioritized_genes
samples = pred_ess.axes[1]

# Configurations
predictions = DataFrame(None, index=genes, columns=samples)
spearman = make_scorer(spearm_cor_func, greater_is_better=True)

X_train_pre = train_exp
X_test_pre = pred_exp

# Filter by coeficient variation
var_thres = VarianceThreshold(0.8).fit(X_train_pre)
X_train_pre = var_thres.transform(X_train_pre)
X_test_pre = var_thres.transform(X_test_pre)

# features_to_keep = X_train_pre.std() / X_train_pre.mean() > 0.1
# X_train_pre = X_train_pre.loc[:, features_to_keep.values]
# X_test_pre = X_test_pre.loc[:, features_to_keep.values]

#X_train_pre = train_exp.join(train_cnv)
#X_test_pre = pred_exp.join(pred_cnv)


for gene in genes:
    # Assemble prediction variables
    X_train = X_train_pre
    y_train = train_ess.ix[:, gene].values
    X_test = X_test_pre

    # Grid search cv
    # pipeline = Pipeline([
    #     ('fs', SelectKBest(f_regression)),
    #     ('clf', PassiveAggressiveRegressor())
    # ])
    #
    # parameters = {
    #     'fs__k': [2000, 2500, 3000, 3500, 4000, 4500],
    #     'clf__epsilon': [1e-12, 1e-9, 1e-4, 1e-3, 1e-2, .1]
    # }
    #
    # clf = GridSearchCV(pipeline, parameters, scoring=spearman)
    # clf.fit(X_train, y_train)
    #
    # best_k = clf.best_estimator_.get_params()['fs__k']
    # best_epsilon = clf.best_estimator_.get_params()['clf__epsilon']
    best_k = 2500
    best_epsilon = 0.01
    best_n_iter = 3

    # Feature selection
    fs = SelectKBest(f_regression, k=best_k).fit(X_train, y_train)
    X_train = fs.transform(X_train_pre)
    X_test = fs.transform(X_test_pre)

    y_preds_test = []
    y_preds_scores = []

    # Training
    cv = ShuffleSplit(len(y_train), n_iter=10)
    for train_i, test_i in cv:
        clf = PassiveAggressiveRegressor(epsilon=best_epsilon, n_iter=best_n_iter).fit(X_train[train_i], y_train[train_i])
        y_preds_scores.append(spearm_cor_func(clf.predict(X_train[test_i]), y_train[test_i]))
        y_preds_test.append(clf.predict(X_test))

    y_preds_scores = Series(y_preds_scores)
    y_preds_test = DataFrame(y_preds_test)

    # Predict
    y_pred = np.mean(y_preds_test[y_preds_scores.notnull()], axis=0).values
    # predictions.loc[gene, :] = y_preds_test.ix[y_preds_scores.argmax()].values

    predictions.loc[gene, :] = y_pred

    if spearmanr(predictions.loc[gene, :], pred_ess.loc[gene])[0] < -0.20:
        print gene

print predictions

# Calculate score
correlations = []
for gene in genes:
    correlations.append(spearmanr(predictions.loc[gene], pred_ess.loc[gene])[0])

# Register run result
score = '%.5f' % np.mean(correlations)
print score