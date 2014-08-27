__author__ = 'emanuel'

import sys
from scipy.stats import spearmanr
from sklearn.linear_model import PassiveAggressiveRegressor
from sklearn.feature_selection import f_regression, SelectKBest
from sklearn.pipeline import Pipeline
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import make_scorer
from dream_2014_functions import read_data_sets


def spearm_cor_func(y_true, y_pred):
    return spearmanr(y_true, y_pred)[0]

# Import data-sets
train_exp, train_cnv, train_ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Setup
spearman = make_scorer(spearm_cor_func, greater_is_better=True)

X_train_pre = train_exp
X_test_pre = leader_exp

# Read gene from output
gene = sys.argv[1]

# Assemble prediction variables
X_train = X_train_pre
y_train = train_ess.ix[:, gene]
X_test = X_test_pre

# Grid search cv
pipeline = Pipeline([
    ('fs', SelectKBest(f_regression)),
    ('clf', PassiveAggressiveRegressor(epsilon=0.01))
])

parameters = {
    'fs__k': [1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250, 3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250]
}

clf = GridSearchCV(pipeline, parameters, scoring=spearman)
clf.fit(X_train, y_train)

best_k = clf.best_estimator_.get_params()['fs__k']

# Feature selection
fs = SelectKBest(f_regression, k=best_k).fit(X_train, y_train)
X_train = fs.transform(X_train)
X_test = fs.transform(X_test)

# Estimation
clf = PassiveAggressiveRegressor(epsilon=0.01)
y_test_pred = clf.fit(X_train, y_train).predict(X_test)

# Store results
output = gene + '\t' + gene + '\t' + '\t'.join([str(x) for x in y_test_pred])

with open('_submissions_temp/' + gene + '.txt', 'w') as f:
    f.write(output)

print gene, best_k