__author__ = 'emanuel'

import numpy as np
import operator
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.metrics import make_scorer
from pandas import DataFrame, Series
from dream_2014_functions import read_data_sets


def spearm_cor_func(expected, pred):
    return spearmanr(expected, pred)[0]

# Import data
exp, cnv, ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Split training data-set in two
train_exp = exp.iloc[range(0, 50), ]
train_cnv = cnv.iloc[range(0, 50), ]
train_ess = ess.iloc[range(0, 50), ]

pred_exp = exp.iloc[range(51, 66), ]
pred_cnv = cnv.iloc[range(51, 66), ]
pred_ess = ess.iloc[range(51, 66), ].T

# Predicted genes
samples = pred_ess.axes[1]

# Configuration
X_train_pre = train_exp
X_test_pre = pred_exp

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
important_features_top_100 = [x for x in important_features.value_counts().head(100).index]

predictions = DataFrame(None, index=prioritized_genes, columns=samples)
my_spearm_cor_func = make_scorer(spearm_cor_func, greater_is_better=True)

# Assemble prediction variables
X_train = X_train_pre.loc[:, important_features_top_100]
X_test = X_test_pre.loc[:, important_features_top_100]

# Normalization
scaler = StandardScaler().fit(X_train)
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)

for gene in prioritized_genes:
    y_train = train_ess.ix[:, gene]

    # Estimation
    clf = RidgeCV(scoring=my_spearm_cor_func)
    y_test_pred = clf.fit(X_train, y_train).predict(X_test)

    print gene, X_train.shape

    # Store results
    predictions.ix[gene] = y_test_pred


# Calculate score
correlations = []
for gene in prioritized_genes:
    correlations.append(spearmanr(predictions.loc[gene], pred_ess.loc[gene])[0])

# Register run result
score = '%.5f' % np.mean(correlations)

print score