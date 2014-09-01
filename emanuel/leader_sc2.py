__author__ = 'emanuel'

import numpy as np
import zipfile
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV, PassiveAggressiveRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_regression
from sklearn.metrics import make_scorer
from sklearn.grid_search import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.cross_validation import ShuffleSplit
from pandas import DataFrame, Series
from dream_2014_functions import read_data_sets, save_gct_data, write_features, submit_solution, ev_code_sc2


def spearm_cor_func(expected, pred):
    return spearmanr(expected, pred)[0]

# Folders
submission_filename_prefix = 'sc2_emanuel_phase2_'

# Import data-sets
train_exp, train_cnv, train_ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Configurations
predictions = DataFrame(None, index=prioritized_genes, columns=leader_exp.axes[0])
spearman = make_scorer(spearm_cor_func, greater_is_better=True)
predictions_features = {}

X_train_pre = train_exp
X_test_pre = leader_exp

# Filter by coeficient variation
var_thres = 0.7
filter_thres = VarianceThreshold(var_thres).fit(X_train_pre)
X_train_pre = X_train_pre.loc[:, filter_thres.get_support()]
X_test_pre = X_test_pre.loc[:, filter_thres.get_support()]

features = X_train_pre.columns.values

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

    print gene, X_train.shape

    # Store results
    predictions.ix[gene] = y_pred
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