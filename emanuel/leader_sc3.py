__author__ = 'emanuel'

import numpy as np
import operator
import zipfile
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV, PassiveAggressiveRegressor
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.cross_validation import ShuffleSplit
from sklearn.metrics import make_scorer
from pandas import DataFrame, Series
from dream_2014_functions import read_data_sets, save_gct_data, submit_solution, write_features_sc3, ev_code_sc3


def spearm_cor_func(expected, pred):
    return spearmanr(expected, pred)[0]

# Folders
submission_filename_prefix = 'sc3_emanuel_phase2_'

# Import data
train_exp, train_cnv, train_ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

X_train_pre = train_exp
X_test_pre = leader_exp

var_thres = VarianceThreshold(0.65).fit(X_train_pre)
X_train_pre = X_train_pre.loc[:, var_thres.get_support()]
X_test_pre = X_test_pre.loc[:, var_thres.get_support()]

# Prepare features
features = X_train_pre.columns
important_features = []

for gene in prioritized_genes:
    # Assemble prediction variables
    X_train = X_train_pre
    y_train = train_ess.ix[:, gene]
    X_test = X_test_pre

    # Feature selection
    fs = SelectKBest(f_regression, k=100)
    X_train = fs.fit_transform(X_train, y_train)
    X_test = fs.transform(X_test)
    gene_features = features[fs.get_support()]

    # Store gene features
    important_features.append(gene_features.values)

    print gene, X_test.shape

# Filter features
important_features = Series([feature for feature_list in important_features for feature in feature_list])
important_features_top_100 = [x for x in important_features.value_counts().head(100).index]

predictions = DataFrame(None, index=prioritized_genes, columns=leader_exp.axes[0])
spearman = make_scorer(spearm_cor_func, greater_is_better=True)

# Assemble prediction variables
X_train = X_train_pre.loc[:, important_features_top_100]
X_test = X_test_pre.loc[:, important_features_top_100]

for gene in prioritized_genes:
    y_train = train_ess.ix[:, gene]

    y_preds_test = []
    y_preds_scores = []

    # Training
    cv = ShuffleSplit(len(y_train), n_iter=5)
    for train_i, test_i in cv:
        clf = PassiveAggressiveRegressor(epsilon=0.01, n_iter=3).fit(X_train.ix[train_i, :], y_train[train_i])
        y_preds_scores.append(spearm_cor_func(clf.predict(X_train.ix[test_i, :]), y_train[test_i]))
        y_preds_test.append(clf.predict(X_test))

    y_preds_scores = Series(y_preds_scores)
    y_preds_test = DataFrame(y_preds_test)

    # Predict
    y_pred = np.mean(y_preds_test[y_preds_scores.notnull()], axis=0).values

    print gene, X_train.shape

    # Store results
    predictions.ix[gene] = y_pred

filename_gct = save_gct_data(predictions, submission_filename_prefix)
print '[DONE]: Saved to file ' + filename_gct

filename_txt = write_features_sc3(important_features_top_100, submission_filename_prefix)
print '[DONE]: Saved to file ' + filename_txt

filename_zip = filename_gct.split('.')[0] + '.zip'

zip_file = zipfile.ZipFile(filename_zip, 'a')
zip_file.write(filename_gct)
zip_file.write(filename_txt)
zip_file.close()

submit_solution(filename_zip, filename_zip.split('/')[1], ev_code_sc3)