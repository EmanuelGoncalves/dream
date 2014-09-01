__author__ = 'emanuel'

import sys
import numpy as np
from pandas import Series, DataFrame
from scipy.stats import spearmanr
from sklearn.linear_model import PassiveAggressiveRegressor
from sklearn.feature_selection import f_regression, SelectKBest
from sklearn.pipeline import Pipeline
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import make_scorer
from sklearn.cross_validation import ShuffleSplit
from dream_2014_functions import read_data_sets


def spearm_cor_func(y_true, y_pred):
    return spearmanr(y_true, y_pred)[0]


# Import data-sets
train_exp, train_cnv, train_ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Setup
spearman = make_scorer(spearm_cor_func, greater_is_better=True)

X_train_pre = train_exp
X_test_pre = leader_exp

# Filter by coeficient variation
features_to_keep = X_train_pre.std() / X_train_pre.mean() > 0.1
X_train_pre = X_train_pre.loc[:, features_to_keep.values]
X_test_pre = X_test_pre.loc[:, features_to_keep.values]

# Read gene from output
gene = sys.argv[1]

# Assemble prediction variables
X_train = X_train_pre
y_train = train_ess.ix[:, gene]
X_test = X_test_pre

# Grid search cv
pipeline = Pipeline([
    ('fs', SelectKBest(f_regression)),
    ('clf', PassiveAggressiveRegressor())
])

parameters = {
    'fs__k': [2000, 2500, 3000, 3500, 4000, 4500],
    'clf__epsilon': [1e-15, 1e-13, 1e-12, 1e-11, 1e-9, 1e-6, 1e-4, 1e-3, 1e-2, .1],
    'clf__n_iter': range(1, 100)
}

clf = GridSearchCV(pipeline, parameters, spearman)
clf.fit(X_train, y_train)

best_k = clf.best_estimator_.get_params()['fs__k']
best_epsilon = clf.best_estimator_.get_params()['clf__epsilon']
best_n_iter = clf.best_estimator_.get_params()['clf__n_iter']

# Feature selection
fs = SelectKBest(f_regression, k=best_k).fit(X_train, y_train)
X_train = fs.transform(X_train)
X_test = fs.transform(X_test)

y_preds_test = []
y_preds_scores = []

cv = ShuffleSplit(len(y_train), n_iter=10)

for train_i, test_i in cv:
    clf = PassiveAggressiveRegressor(epsilon=best_epsilon, n_iter=best_n_iter)

    clf.fit(X_train[train_i], y_train[train_i])

    y_preds_scores.append(spearm_cor_func(clf.predict(X_train[test_i]), y_train[test_i]))

    y_preds_test.append(clf.predict(X_test))

y_preds_scores = Series(y_preds_scores)
y_preds_test = DataFrame(y_preds_test)

# Predict
y_test_pred = np.mean(y_preds_test[y_preds_scores.notnull()], axis=0).values

# Store results
output = gene + '\t' + gene + '\t' + '\t'.join([str(x) for x in y_test_pred])

print '[INFO]: ', gene, best_k, best_epsilon, best_n_iter

with open('_submissions_temp/' + gene + '.txt', 'w') as f:
    f.write(output)