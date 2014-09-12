__author__ = 'emanuel'

import numpy as np
from pandas import DataFrame, Series
from scipy.stats import spearmanr
from sklearn.linear_model import PassiveAggressiveRegressor
from sklearn.feature_selection import f_regression, SelectKBest
from sklearn.feature_selection import VarianceThreshold
from sklearn.cross_validation import ShuffleSplit
from sklearn.metrics import make_scorer
from dream_2014_functions import read_data_sets, save_gct_data, submit_solution, ev_code_sc1


def spearm_cor_func(y_true, y_pred):
    return spearmanr(y_true, y_pred)[0]

best_epsilon = 0.01
best_n_iter = 3
best_k = 3500
best_var = 0.65

# Folders
submission_filename_prefix = 'sc1_emanuel_phase2_'

# Import data-sets
train_exp, train_cnv, train_ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Setup
genes = train_ess.axes[1]
samples = leader_exp.axes[0]
predictions = DataFrame(None, index=genes, columns=samples)
spearman = make_scorer(spearm_cor_func, greater_is_better=True)

X_train_pre = train_exp
X_test_pre = leader_exp

# Filter by coeficient variation
var_thres = VarianceThreshold(best_var).fit(X_train_pre)
X_train_pre = var_thres.transform(X_train_pre)
X_test_pre = var_thres.transform(X_test_pre)

for gene in genes:
    # Assemble prediction variables
    X_train = X_train_pre
    y_train = train_ess.ix[:, gene]
    X_test = X_test_pre

    # Feature selection
    fs = SelectKBest(f_regression, k=best_k).fit(X_train, y_train)
    X_train = fs.transform(X_train)
    X_test = fs.transform(X_test)

    # Estimation
    clf = PassiveAggressiveRegressor(epsilon=best_epsilon, n_iter=best_n_iter).fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    # Store results
    predictions.ix[gene] = y_pred

    print gene

filename = save_gct_data(predictions, submission_filename_prefix)
print '[DONE]: Saved to file ' + filename

submit_solution(filename, filename.split('/')[1], ev_code_sc1)
print '[SUBMITED]'