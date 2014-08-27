__author__ = 'emanuel'

from pandas import DataFrame
from scipy.stats import spearmanr
from sklearn.linear_model import PassiveAggressiveRegressor
from sklearn.feature_selection import f_regression, SelectKBest
from sklearn.pipeline import Pipeline
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import make_scorer
from dream_2014_functions import read_data_sets, save_gct_data, submit_solution, ev_code_sc1


def spearm_cor_func(y_true, y_pred):
    return spearmanr(y_true, y_pred)[0]

# Folders
submission_filename_prefix = 'sc1_emanuel_phase2_'

# Import data-sets
train_exp, train_cnv, train_ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

# Setup
genes = train_ess.axes[1]
samples = leader_exp.axes[0]
predictions = DataFrame(None, index=genes, columns=samples)
spearman = make_scorer(spearm_cor_func, greater_is_better=True)

X_train_pre = train_exp
X_test_pre = leader_exp

for gene in genes:
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
        'fs__k': [1500, 2000, 2500, 3000, 3500, 4000]
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
    predictions.ix[gene] = y_test_pred

    print gene, best_k

filename = save_gct_data(predictions, submission_filename_prefix)
print '[DONE]: Saved to file ' + filename

submit_solution(filename, filename.split('/')[1], ev_code_sc1)
print '[SUBMITED]'