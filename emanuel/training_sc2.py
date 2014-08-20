__author__ = 'emanuel'

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_regression
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import make_scorer
from pandas import DataFrame
from dream_2014_functions import read_data_sets


def spearm_cor_func(expected, pred):
    return spearmanr(expected, pred)[0]

# Import data-sets
exp, cnv, ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Split training data-set in two
train_exp = exp.iloc[range(0, 50), ]
train_cnv = cnv.iloc[range(0, 50), ]
train_ess = ess.iloc[range(0, 50), ]

pred_exp = exp.iloc[range(51, 66), ]
pred_cnv = cnv.iloc[range(51, 66), ]
pred_ess = ess.iloc[range(51, 66), ].T

X_train_pre = train_exp
X_test_pre = pred_exp

samples = pred_ess.axes[1]
features = X_train_pre.axes[1]

# Configurations
predictions = DataFrame(None, index=prioritized_genes, columns=samples)

my_spearm_cor_func = make_scorer(spearm_cor_func, greater_is_better=True)

# Filter by coeficient variation
features_to_keep = X_train_pre.std() / X_train_pre.mean() > 0.1
X_train_pre = X_train_pre.loc[:, features_to_keep.values]
X_test_pre = X_test_pre.loc[:, features_to_keep.values]

for gene in prioritized_genes:
    # Assemble prediction variables
    X_train = X_train_pre
    y_train = train_ess.ix[:, gene]
    X_test = X_test_pre

    # Normalization
    scaler = StandardScaler().fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)

    # Feature selection
    clf = RidgeCV()

    parameters = {
        'alphas': [
            [.0001, 0.001, .01, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1],
            [.2, .3, .4, .5, .6],
            [.0001, 0.001, .01, .1],
            [.5, .6, .7, .8, .9, 1],
            [.65, .7, .75, .8, .85, .9, .95, .99, 1],
            [.5, .6, .7, .75, .8, .85, .9, .95, .99, 1]
        ]
    }

    gs_clf = GridSearchCV(clf, parameters)
    gs_clf = gs_clf.fit(X_train, y_train)

    best_parameters, score, _ = max(gs_clf.grid_scores_, key=lambda x: x[1])
    print best_parameters, score

    # Estimation
    clf = RidgeCV(alphas=best_parameters['alphas'], gcv_mode='auto')
    y_test_pred = clf.fit(X_train, y_train).predict(X_test)

    print gene, X_train.shape

    # Store results
    predictions.ix[gene] = y_test_pred

# Calculate score
correlations = []
for gene in prioritized_genes:
    correlations.append(spearmanr(predictions.loc[gene], pred_ess.loc[gene])[0])

# Register run result
score = '%.5f' % np.mean(correlations)

print score