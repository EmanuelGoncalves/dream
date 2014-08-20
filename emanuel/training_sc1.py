__author__ = 'emanuel'

# Set-up workspace
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV, PassiveAggressiveRegressor
from sklearn.metrics import make_scorer
from pandas import DataFrame, Series
from dream_2014_functions import read_data_sets, read_annotations


def hill_function(matrix, hill_coef=6):
    return 1 / ((np.median(exp, axis=0) / matrix) ** hill_coef + 1)


def spearm_cor_func(expected, pred):
    return spearmanr(expected, pred)[0]

# Import data-sets
exp, cnv, ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Import annotations
gene_annot = read_annotations()

# Preprocess expression data-sets
exp = exp.drop(gene_annot[gene_annot.isnull().values].axes[0].values, 1)
exp = hill_function(exp)

# Split training data-set in two
train_exp = exp.iloc[range(50), ]
train_cnv = cnv.iloc[range(50), ]
train_ess = ess.iloc[range(50), ]

pred_exp = exp.iloc[range(51, 66), ]
pred_cnv = cnv.iloc[range(51, 66), ]
pred_ess = ess.iloc[range(51, 66), ].T

# Predicted genes
genes = pred_ess.axes[0]
#genes = prioritized_genes
samples = pred_ess.axes[1]

# Configurations
predictions = DataFrame(None, index=genes, columns=samples)
my_spearm_cor_func = make_scorer(spearm_cor_func, greater_is_better=True)
cv_thres = 0.3

X_train_pre = train_exp
X_test_pre = pred_exp

# Filter by coeficient variation
features_to_keep = X_train_pre.std() / X_train_pre.mean() > cv_thres
X_train_pre = X_train_pre.loc[:, features_to_keep.values]
X_test_pre = X_test_pre.loc[:, features_to_keep.values]

for gene in genes:
    # Assemble prediction variables
    X_train = X_train_pre
    y_train = train_ess.ix[:, gene]
    X_test = X_test_pre

    # Estimation
    clf = PassiveAggressiveRegressor(epsilon=0.01)
    y_test_pred = clf.fit(X_train, y_train).predict(X_test)

    # Store results
    predictions.ix[gene] = y_test_pred

    print gene, X_train.shape

# Calculate score
correlations = []
for gene in genes:
    correlations.append(spearmanr(predictions.loc[gene], pred_ess.loc[gene])[0])

# Register run result
score = '%.5f' % np.mean(correlations)
config_output = 'CV+HillFunction\t'+str(cv_thres)+',6\tPassiveAggressiveRegressor+Clustering\t\t' + score
print config_output

with open('emanuel/training_sc1.log', 'a') as f:
    f.write(config_output + '\n')

