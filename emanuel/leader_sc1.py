__author__ = 'emanuel'

import numpy as np
from sklearn.linear_model import PassiveAggressiveRegressor
from pandas import DataFrame
from dream_2014_functions import read_data_sets, save_gct_data, submit_solution, ev_code_sc1, read_annotations


def hill_function(matrix, hill_coef=9):
    return 1 / ((np.median(matrix, axis=0) / matrix) ** hill_coef + 1)

# Folders
submission_filename_prefix = 'sc1_emanuel_phase2_'

# Import data-sets
train_exp, train_cnv, train_ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Import annotations
gene_annot = read_annotations()

# Preprocess expression data-sets
probes_remove = gene_annot[gene_annot.isnull().values].axes[0].values

train_exp = train_exp.drop(probes_remove, 1)
train_exp = hill_function(train_exp)

leader_exp = leader_exp.drop(probes_remove, 1)
leader_exp = hill_function(leader_exp)

# Setup
genes = train_ess.axes[1]
samples = leader_exp.axes[0]
predictions = DataFrame(None, index=genes, columns=samples)
cv_thres = 0.25

X_train_pre = train_exp
X_test_pre = leader_exp

# Filter by coeficient variation
features_to_keep = X_train_pre.std() / X_train_pre.mean() > cv_thres
X_train_pre = X_train_pre.loc[:, features_to_keep.values]
X_test_pre = X_test_pre.loc[:, features_to_keep.values]

for gene in genes:
    # Assemble prediction variables
    X_train = X_train_pre
    y_train = train_ess.ix[:, gene]
    X_test = X_test_pre

    # Estimation
    clf = PassiveAggressiveRegressor(epsilon=0.01)
    y_test_pred = clf.fit(X_train, y_train).predict(X_test)

    # Store results
    predictions.ix[gene] = y_test_pred

    print gene, X_train.shape

print '[INFO]: Coeficient Variation threshold: ' + str(cv_thres)

filename = save_gct_data(predictions, submission_filename_prefix)
print '[DONE]: Saved to file ' + filename

submit_solution(filename, filename.split('/')[1], ev_code_sc1)
print '[SUBMITED]'