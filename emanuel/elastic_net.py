__author__ = 'emanuel'

from pandas import read_csv, DataFrame
from sklearn.linear_model import ElasticNet, PassiveAggressiveRegressor, LinearRegression
from sklearn.preprocessing import StandardScaler
import numpy as np

# Data-sets files
train_exp_file = 'data/CCLE_expression_training.gct'
train_cnv_file = 'data/CCLE_copynumber_training.gct'
train_ess_file = 'data/Achilles_v2.9_training.gct'

leader_exp_file = 'data/CCLE_expression_leaderboard.gct'
leader_cnv_file = 'data/CCLE_copynumber_leaderboard.gct'

# Folders
submissions_folder = 'submissions/'
submission_filename = 'umebi_emanuel_20.gct'

# Utilities function
def read_gct(filename):
    data = read_csv(filename, sep='\t', header=2, index_col=0)
    del data['Description']
    return data.T

def save_gct_data(dataframe, filename=submission_filename, folder=submissions_folder, sep='\t'):
    rows = dataframe.shape[0]
    cols = dataframe.shape[1]

    with open(folder + filename, 'w') as f:
        f.write('#1.2\n')
        f.write(str(rows) + '\t' + str(cols) + '\n')

        f.write('Name' + sep + 'Description')
        for header in dataframe.axes[1]:
            f.write(sep + header)
        f.write('\n')

        for i in range(rows):
            row_name = dataframe.axes[0][i]
            f.write(row_name + sep + row_name)

            for j in range(cols):
                f.write(sep + '%.5f' % dataframe.ix[i, j])

            f.write('\n')

# Methods
methods = ['ElasticNet', 'PassiveAggressiveRegressor', 'LinearRegression']

# Configurations
method = methods[2]
filter_thres = 95

print 'Method: ' + str(method)
print 'Threshold: ' + str(filter_thres)

# Import data
train_exp = read_gct(train_exp_file)
train_cnv = read_gct(train_cnv_file)
train_ess = read_gct(train_ess_file)

leader_exp = read_gct(leader_exp_file)
leader_cnv = read_gct(leader_cnv_file)

genes = train_ess.axes[1]

predictions = DataFrame(None, index=genes, columns=leader_exp.axes[0])

# Assemble trainning and prediction features datasets
train_features = train_exp.join(train_cnv)
leader_features = leader_exp.join(leader_cnv)

# Filter features by variance
train_features_var = np.var(train_features.ix[:, ])
train_features_var_thres = np.percentile(train_features_var, filter_thres)

train_features = train_features.ix[:, train_features_var > train_features_var_thres]
leader_features = leader_features.ix[:, train_features_var > train_features_var_thres]

# Standardize features
scaler = StandardScaler().fit(train_features)
train_features = scaler.transform(train_features)
leader_features = scaler.transform(leader_features)

# Perform elastic net
for gene in genes:
    print gene

    # Assemble prediction variables
    X_train = train_features
    y_train = train_ess.ix[:, gene]
    X_test = leader_features

    # Perform method
    if method == methods[0]:
        en = ElasticNet(alpha=0.1, l1_ratio=0.7, warm_start=True)
        y_test_pred = en.fit(X_train, y_train).predict(X_test)

    elif method == methods[1]:
        par = PassiveAggressiveRegressor()
        y_test_pred = par.fit(X_train, y_train).predict(X_test)

    elif method == methods[2]:
        lr = LinearRegression()
        y_test_pred = lr.fit(X_train, y_train).predict(X_test)

    # Store results
    predictions.ix[gene] = y_test_pred

save_gct_data(predictions)

print '[DONE]: Saved to file ' + submission_filename