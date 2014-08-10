__author__ = 'emanuel'

# Set-up workspace
import sys
sys.path.append('/Users/emanuel/Documents/projects_data_analysis/dream/emanuel/')

from pandas import DataFrame
from sklearn.linear_model import PassiveAggressiveRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from dream_2014_functions import read_gct, save_gct_data, submit_solution

# Data-sets files
train_exp_file = 'data/CCLE_expression_training.gct'
train_cnv_file = 'data/CCLE_copynumber_training.gct'
train_ess_file = 'data/Achilles_v2.9_training.gct'

leader_exp_file = 'data/CCLE_expression_leaderboard.gct'
leader_cnv_file = 'data/CCLE_copynumber_leaderboard.gct'

# Folders
submissions_folder = 'submissions/'
submission_filename_prefix = 'sc1_emanuel_'

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

genes = train_ess.axes[1]

predictions = DataFrame(None, index=genes, columns=leader_exp.axes[0])

# Assemble trainning and prediction features datasets
train_features = train_exp
leader_features = leader_exp

# Assemble predictions
X_train_pre = train_features
X_test_pre = leader_features

var_fs = VarianceThreshold(0.25)
X_train_pre = var_fs.fit_transform(X_train_pre)
X_test_pre = var_fs.transform(X_test_pre)

# Perform elastic net
for gene in genes:
    # Assemble prediction variables
    y_train = train_ess.ix[:, gene]
    X_train = X_train_pre
    X_test = X_test_pre

    # Standardize features
    X_train = StandardScaler().fit_transform(X_train)
    X_test = StandardScaler().fit_transform(X_test)

    print gene, X_train.shape

    # Training
    par = PassiveAggressiveRegressor(epsilon=0.01)
    par.fit(X_train, y_train)
    y_test_pred = par.predict(X_test)

    # Store results
    predictions.ix[gene] = y_test_pred

filename = save_gct_data(predictions, submission_filename_prefix)
print '[DONE]: Saved to file ' + filename

submit_solution(filename, filename.split('/')[1], ev_code_sc1)
print '[SUBMITED]'