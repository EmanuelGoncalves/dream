__author__ = 'emanuel'

import numpy as np
import zipfile
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV, PassiveAggressiveRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_regression
from sklearn.metrics import make_scorer
from sklearn.grid_search import GridSearchCV
from pandas import DataFrame
from dream_2014_functions import read_data_sets, save_gct_data, write_features, submit_solution, ev_code_sc2


def spearm_cor_func(expected, pred):
    return spearmanr(expected, pred)[0]

# Folders
submission_filename_prefix = 'sc2_emanuel_phase2_'

# Import data-sets
train_exp, train_cnv, train_ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Configurations
predictions = DataFrame(None, index=prioritized_genes, columns=leader_exp.axes[0])
predictions_features = {}

X_train_pre = train_exp
X_test_pre = leader_exp

features = X_train_pre.axes[1]

for gene in prioritized_genes:
    # Assemble prediction variables
    X_train = X_train_pre
    y_train = train_ess.ix[:, gene]
    X_test = X_test_pre

    # Feature selection
    fs = SelectKBest(f_regression).fit(X_train, y_train)
    X_train = fs.transform(X_train)
    X_test = fs.transform(X_test)
    gene_features = features[fs.get_support()]

    # Estimation
    #clf = RidgeCV(gcv_mode='auto')
    clf = PassiveAggressiveRegressor(epsilon=0.01)
    y_test_pred = clf.fit(X_train, y_train).predict(X_test)

    print gene, X_train.shape

    # Store results
    predictions.ix[gene] = y_test_pred
    predictions_features[gene] = gene_features

filename_gct = save_gct_data(predictions, submission_filename_prefix)
print '[DONE]: Saved to file ' + filename_gct

filename_txt = write_features(predictions_features, submission_filename_prefix)
print '[DONE]: Saved to file ' + filename_txt

filename_zip = filename_gct.split('.')[0] + '.zip'

zip_file = zipfile.ZipFile(filename_zip, 'a')
zip_file.write(filename_gct)
zip_file.write(filename_txt)
zip_file.close()

submit_solution(filename_zip, filename_zip.split('/')[1], ev_code_sc2)