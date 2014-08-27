__author__ = 'emanuel'

import numpy as np
import zipfile
from sklearn.linear_model import PassiveAggressiveRegressor, RidgeCV
from sklearn.feature_selection import SelectKBest, f_regression
from pandas import DataFrame, Series
from dream_2014_functions import read_data_sets, save_gct_data, submit_solution, write_features_sc3, ev_code_sc3

# Folders
submission_filename_prefix = 'sc3_emanuel_phase2_'

# Import data
train_exp, train_cnv, train_ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

X_train_pre = train_exp
X_test_pre = leader_exp

# Prepare features
features = X_train_pre.axes[1]
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
important_features_top_100 = important_features.value_counts().head(100)

predictions = DataFrame(None, index=prioritized_genes, columns=leader_exp.axes[0])

for gene in prioritized_genes:
    # Assemble prediction variables
    X_train = X_train_pre.loc[:, important_features_top_100.index]
    y_train = train_ess.ix[:, gene]
    X_test = X_test_pre.loc[:, important_features_top_100.index]

    # Estimation
    #clf = PassiveAggressiveRegressor(epsilon=0.01)
    clf = RidgeCV()
    y_test_pred = clf.fit(X_train, y_train).predict(X_test)

    print gene, X_train.shape

    # Store results
    predictions.ix[gene] = y_test_pred

filename_gct = save_gct_data(predictions, submission_filename_prefix)
print '[DONE]: Saved to file ' + filename_gct

filename_txt = write_features_sc3(important_features_top_100.index, submission_filename_prefix)
print '[DONE]: Saved to file ' + filename_txt

filename_zip = filename_gct.split('.')[0] + '.zip'

zip_file = zipfile.ZipFile(filename_zip, 'a')
zip_file.write(filename_gct)
zip_file.write(filename_txt)
zip_file.close()

submit_solution(filename_zip, filename_zip.split('/')[1], ev_code_sc3)