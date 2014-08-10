__author__ = 'emanuel'

# Set-up workspace
import sys
sys.path.append('/Users/emanuel/Documents/projects_data_analysis/dream/emanuel/')

import numpy as np
import zipfile
import operator
from sklearn.linear_model import ElasticNetCV, RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_regression
from pandas import DataFrame
from dream_2014_functions import read_gct, save_gct_data, read_features, submit_solution, write_features_sc3

# Data-sets files
train_exp_file = 'data/CCLE_expression_training.gct'
train_cnv_file = 'data/CCLE_copynumber_training.gct'
train_ess_file = 'data/Achilles_v2.9_training.gct'

leader_exp_file = 'data/CCLE_expression_leaderboard.gct'
leader_cnv_file = 'data/CCLE_copynumber_leaderboard.gct'

gene_list = 'data/prioritized_gene_list.txt'

# Folders
submissions_folder = 'submissions/'
submission_filename_prefix = 'sc3_emanuel_'

# Evaluation codes
ev_code_sc1 = 2468319
ev_code_sc2 = 2468322
ev_code_sc3 = 2482339

# Import data
train_exp = read_gct(train_exp_file)
train_cnv = read_gct(train_cnv_file)
train_ess = read_gct(train_ess_file)

leader_exp = read_gct(leader_exp_file)
leader_cnv = read_gct(leader_cnv_file)

# Predicted genes
genes = np.genfromtxt(gene_list, dtype='str')
samples = leader_exp.axes[0]

# Configurations
predictions = DataFrame(None, index=genes, columns=samples)

var_fs_thres = 0.25

X_train_pre = train_exp.join(train_cnv)
X_test_pre = leader_exp.join(leader_cnv)
features = X_train_pre.axes[1]

# Feature selection
features_dict = read_features('submissions/sc2_emanuel_3.txt')
features_soreted = sorted(features_dict.iteritems(), key=operator.itemgetter(1), reverse=True)
features = [features_soreted[i][0] for i in range(100)]

for gene in genes:
    # Assemble prediction variables
    X_train = X_train_pre
    y_train = train_ess.ix[:, gene]
    X_test = X_test_pre

    # Normalization
    scaler = StandardScaler().fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)

    # Estimation
    #clf = ElasticNetCV(l1_ratio=[.1, .12, .15, .17, .2, .25, .3, .5, .7, .9, 1])
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