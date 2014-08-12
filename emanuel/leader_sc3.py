__author__ = 'emanuel'

import numpy as np
import zipfile
import operator
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from pandas import DataFrame
from dream_2014_functions import read_data_sets, save_gct_data, read_features, submit_solution, write_features_sc3, ev_code_sc3

# Folders
submission_filename_prefix = 'sc3_emanuel_phase2_'

# Import data
train_exp, train_cnv, train_ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Predicted genes
samples = leader_exp.axes[0]

# Configurations
predictions = DataFrame(None, index=prioritized_genes, columns=samples)
var_fs_thres = 0.25
cv_n = 5

# Feature selection
features_dict = read_features('submissions/sc2_emanuel_3.txt')
features_soreted = sorted(features_dict.iteritems(), key=operator.itemgetter(1), reverse=True)
features = [features_soreted[i][0] for i in range(100)]

X_train_pre = train_exp[features]
X_test_pre = leader_exp[features]

for gene in prioritized_genes:
    # Assemble prediction variables
    X_train = X_train_pre
    y_train = train_ess.ix[:, gene]
    X_test = X_test_pre

    # Normalization
    scaler = StandardScaler().fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)

    # Estimation
    clf = RidgeCV(gcv_mode='auto')
    y_test_pred = clf.fit(X_train, y_train).predict(X_test)

    print gene, X_train.shape, clf.alpha_

    # Store results
    predictions.ix[gene] = y_test_pred

filename_gct = save_gct_data(predictions, submission_filename_prefix)
print '[DONE]: Saved to file ' + filename_gct

filename_txt = write_features_sc3(features, submission_filename_prefix)
print '[DONE]: Saved to file ' + filename_txt

filename_zip = filename_gct.split('.')[0] + '.zip'

zip_file = zipfile.ZipFile(filename_zip, 'a')
zip_file.write(filename_gct)
zip_file.write(filename_txt)
zip_file.close()

submit_solution(filename_zip, filename_zip.split('/')[1], ev_code_sc3)