__author__ = 'emanuel'

# Set-up workspace
import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFECV
from sklearn.metrics import make_scorer
from pandas import DataFrame
from dream_2014_functions import read_data_sets


def spearm_cor_func(expected, pred):
    return spearmanr(expected, pred)[0]

# Import data-sets
exp, cnv, ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Split training data-set in two
train_exp = exp.iloc[range(50), ]
train_cnv = cnv.iloc[range(50), ]
train_ess = ess.iloc[range(50), ]

pred_exp = exp.iloc[range(51, 66), ]
pred_cnv = cnv.iloc[range(51, 66), ]
pred_ess = ess.iloc[range(51, 66), ].T

# Predicted genes
genes = pred_ess.axes[0]
samples = pred_ess.axes[1]

# Configurations
predictions = DataFrame(None, index=genes, columns=samples)
var_fs_thres = 0.18
par_epsilon = 0.01

my_spearm_cor_func = make_scorer(spearm_cor_func, greater_is_better=True)

X_train_pre = train_exp
X_test_pre = pred_exp

for gene in genes:
    # Assemble prediction variables
    X_train = X_train_pre
    y_train = train_ess.ix[:, gene]
    X_test = X_test_pre

    # Normalization
    scaler = StandardScaler().fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)

    # Feature selection
    clf = RidgeCV()

    fs = RFECV(clf, step=150, scoring=my_spearm_cor_func)
    X_train = fs.fit_transform(X_train, y_train)
    X_test = fs.transform(X_test)

    print gene, X_train.shape

    # Estimation
    y_test_pred = clf.fit(X_train, y_train).predict(X_test)

    # Store results
    predictions.ix[gene] = y_test_pred

# Calculate score
correlations = []
for gene in genes:
    correlations.append(spearmanr(predictions.loc[gene], pred_ess.loc[gene])[0])

# Register run result
score = '%.5f' % np.mean(correlations)
config_output = 'RFECV\t300\tspearman\tRidgeCV\t\t' + score
print config_output

with open('emanuel/training_sc1.log', 'a') as f:
    f.write(config_output + '\n')

