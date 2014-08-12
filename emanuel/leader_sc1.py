__author__ = 'emanuel'

# Set-up workspace
import sys
sys.path.append('/Users/emanuel/Documents/projects_data_analysis/dream/emanuel/')

from pandas import DataFrame
from sklearn.linear_model import PassiveAggressiveRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_regression
from dream_2014_functions import read_data_sets, save_gct_data, submit_solution, ev_code_sc1

# Folders
submission_filename_prefix = 'sc1_emanuel_phase2_'

# Import data-sets
train_exp, train_cnv, train_ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

genes = train_ess.axes[1]
samples = leader_exp.axes[0]

predictions = DataFrame(None, index=genes, columns=samples)

# Assemble trainning and prediction features datasets
train_features = train_exp
leader_features = leader_exp

# Assemble predictions
X_train_pre = train_features
X_test_pre = leader_features

# Perform elastic net
for gene in genes:
    # Assemble prediction variables
    y_train = train_ess.ix[:, gene]
    X_train = X_train_pre
    X_test = X_test_pre

    # Standardize features
    X_train = StandardScaler().fit_transform(X_train)
    X_test = StandardScaler().fit_transform(X_test)

    # Feature selection
    fs = SelectKBest(f_regression, k=5000)
    X_train = fs.fit_transform(X_train, y_train)
    X_test = fs.transform(X_test)

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