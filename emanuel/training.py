__author__ = 'emanuel'

# Set-up workspace
import sys
sys.path.append('/Users/emanuel/Documents/projects_data_analysis/dream/emanuel/')

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import ElasticNet, PassiveAggressiveRegressor, LinearRegression, ElasticNetCV, enet_path
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectPercentile, f_regression, SelectKBest, RFE, VarianceThreshold, GenericUnivariateSelect
from sklearn.ensemble import ExtraTreesClassifier, GradientBoostingRegressor
from sklearn import svm
from pandas import DataFrame
from dream_2014_functions import read_gct


def register_trainning(method, features='NA', normalization='NA', feature_selection='NA', feature_selection_thres='NA', score='NA', others='NA', sep='\t'):
    with open('emanuel/training.log', 'a') as f:
        f.write(method + sep + features + sep + normalization + sep + feature_selection + sep + feature_selection_thres + sep + score + sep + others + '\n')

# Data-sets files
train_exp_file = 'data/CCLE_expression_training.gct'
train_cnv_file = 'data/CCLE_copynumber_training.gct'
train_ess_file = 'data/Achilles_v2.9_training.gct'

# Methods
methods = ['ElasticNet', 'PassiveAggressiveRegressor', 'LinearRegression', 'SVM']
feature_selection_methods = ['Percentil', 'ExtraTrees', 'SelectKBest', 'Variance', 'Variance2', 'GradientBoostingRegressor']

# Import data
exp = read_gct(train_exp_file)
cnv = read_gct(train_cnv_file)
ess = read_gct(train_ess_file)

# Split training data-set in two
train_exp = exp.iloc[range(0, 25), ]
train_cnv = cnv.iloc[range(0, 25), ]
train_ess = ess.iloc[range(0, 25), ]

pred_exp = exp.iloc[range(26, 45), ]
pred_cnv = cnv.iloc[range(26, 45), ]
pred_ess = ess.iloc[range(26, 45), ].T

# Predicted genes
genes = pred_ess.axes[0]
samples = pred_ess.axes[1]

# Configurations
predictions = DataFrame(None, index=genes, columns=samples)

method = methods[1]
feature_selection_method = 'NA'

feat_sel_thres = 1500

alpha = 0.25
l1_ratio = 0.9
others = 'NA'

X_train = train_exp
X_test = pred_exp

fs = GradientBoostingRegressor(loss='huber', n_estimators=feat_sel_thres)
X_train = fs.fit_transform(X_train, train_ess.ix[:, genes[4527]])
X_test = fs.transform(X_test)

others = 'huber; '+str(feat_sel_thres)

for gene in genes:
    # Assemble prediction variables
    y_train = train_ess.ix[:, gene]

    # Feature selection
    if feature_selection_method == feature_selection_methods[0]:
        fs = SelectPercentile(f_regression, feat_sel_thres).fit(X_train, y_train)
        X_train = fs.transform(X_train)
        X_test = fs.transform(X_test)

    elif feature_selection_method == feature_selection_methods[1]:
        fs = ExtraTreesClassifier(n_estimators=feat_sel_thres, bootstrap=False)
        X_train = fs.fit_transform(X_train, y_train)
        X_test = fs.transform(X_test)
        others = 'max_features=0.8; bootstrap=True'

    elif feature_selection_method == feature_selection_methods[2]:
        k_best = int(train_exp.shape[1] * feat_sel_thres)
        fs = SelectKBest(f_regression, k=k_best).fit(X_train, y_train)
        X_train = fs.transform(X_train)
        X_test = fs.transform(X_test)

    elif feature_selection_method == feature_selection_methods[3]:
        train_features_var = np.var(X_train.ix[:, ])
        train_features_var_thres = np.percentile(train_features_var, feat_sel_thres)
        X_train = X_train.ix[:, train_features_var > train_features_var_thres]
        X_test = X_test.ix[:, train_features_var > train_features_var_thres]

    elif feature_selection_method == feature_selection_methods[4]:
        fs_train = VarianceThreshold(threshold=feat_sel_thres).fit(X_train)
        X_train = fs_train.transform(X_train)
        X_test = fs_train.transform(X_test)

        fs_test = VarianceThreshold(threshold=feat_sel_thres).fit(X_test)
        X_train = fs_test.transform(X_train)
        X_test = fs_test.transform(X_test)

        others = '2 Variance filtering'

    elif feature_selection_method == feature_selection_methods[5]:
        fs = GradientBoostingRegressor(loss='lad', n_estimators=feat_sel_thres)
        X_train = fs.fit_transform(X_train, y_train)
        X_test = fs.transform(X_test)

    # Normalization
    scaler = StandardScaler().fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)

    print gene, X_train.shape

    # Perform method
    if method == methods[0]:
        en = ElasticNet(alpha=alpha, l1_ratio=l1_ratio)
        y_test_pred = en.fit(X_train, y_train).predict(X_test)

        others = 'alpha='+str(alpha)+';'+'l1_ratio'+str(l1_ratio)

    elif method == methods[1]:
        par = PassiveAggressiveRegressor()
        y_test_pred = par.fit(X_train, y_train).predict(X_test)

    elif method == methods[2]:
        lr = LinearRegression()
        y_test_pred = lr.fit(X_train, y_train).predict(X_test)

    elif method == methods[3]:
        clf = svm.SVR()
        clf.fit(X_train, y_train)
        y_test_pred = clf.predict(X_test)

    # Store results
    predictions.ix[gene] = y_test_pred

# Calculate score
correlations = []
for gene in genes:
    correlations.append(spearmanr(predictions.loc[gene], pred_ess.loc[gene])[0])

# Register run result
score = '%.5f' % np.mean(correlations)

print score
register_trainning(
    method,
    features='GEX',
    normalization='u=0;s=1',
    feature_selection=feature_selection_method,
    feature_selection_thres=str(feat_sel_thres),
    score=score,
    others=others
)


