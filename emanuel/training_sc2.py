__author__ = 'emanuel'

import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from sklearn.linear_model import *
from sklearn.svm import *
from sklearn.pipeline import *
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_regression
from sklearn.grid_search import GridSearchCV
from sklearn.cross_validation import ShuffleSplit
from sklearn.ensemble import *
from sklearn.tree import *
from sklearn.metrics import make_scorer
from pandas import DataFrame, Series
from dream_2014_functions import read_data_sets


def spearm_cor_func(expected, pred):
    return spearmanr(expected, pred)[0]


def hill_function(matrix, hill_coef=2):
    return 1 / ((np.median(matrix, axis=0) / matrix) ** hill_coef + 1)

# Import data-sets
exp, cnv, ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Split training data-set in two
train_exp = exp.iloc[range(26, 66), ]
train_cnv = cnv.iloc[range(26, 66), ]
train_ess = ess.iloc[range(26, 66), ]

pred_exp = exp.iloc[range(25), ]
pred_cnv = cnv.iloc[range(25), ]
pred_ess = ess.iloc[range(25), ].T

X_train_pre = train_exp
X_test_pre = pred_exp

# Configurations
predictions = DataFrame(None, index=prioritized_genes, columns=pred_ess.axes[1])
spearman = make_scorer(spearm_cor_func, greater_is_better=True)

# Filter by coeficient variation
var_thres = 0.65
filter_thres = VarianceThreshold(var_thres).fit(X_train_pre)
X_train_pre = filter_thres.transform(X_train_pre)
X_test_pre = filter_thres.transform(X_test_pre)

for gene in prioritized_genes:
    # Assemble prediction variables
    X_train = X_train_pre
    y_train = train_ess.ix[:, gene].values
    X_test = X_test_pre

    fs = SelectKBest(f_regression).fit(X_train, y_train)
    X_train = fs.transform(X_train)
    X_test = fs.transform(X_test)

    y_preds_test = []
    y_preds_scores = []

    # Training
    cv = ShuffleSplit(len(y_train), n_iter=5)
    for train_i, test_i in cv:
        clf = RidgeCV(gcv_mode='auto').fit(X_train[train_i], y_train[train_i])
        y_preds_scores.append(spearm_cor_func(clf.predict(X_train[test_i]), y_train[test_i]))
        y_preds_test.append(clf.predict(X_test))

    y_preds_scores = Series(y_preds_scores)
    y_preds_test = DataFrame(y_preds_test)

    # Predict
    y_pred = np.mean(y_preds_test[y_preds_scores.notnull()], axis=0).values

    predictions.ix[gene] = y_pred

    print gene, X_train.shape, spearmanr(predictions.loc[gene, :], pred_ess.loc[gene])

# Calculate score
correlations = []
for gene in prioritized_genes:
    correlations.append(spearmanr(predictions.loc[gene], pred_ess.loc[gene])[0])

# Register run result
score = '%.5f' % np.mean(correlations)

print score