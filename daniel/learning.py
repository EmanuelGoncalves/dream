__author__ = 'daniel'

from pandas import DataFrame
from numpy import mean
from scipy.stats import spearmanr
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet, LassoCV, ElasticNetCV, \
    BayesianRidge, PassiveAggressiveRegressor, ARDRegression, LogisticRegression
from sklearn.gaussian_process import GaussianProcess
from sklearn.cross_decomposition import PLSRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.isotonic import IsotonicRegression
from sklearn.feature_selection import RFE
from time import time

from datasets import load_datasets, load_cell_lines, save_gct_data, submit_to_challenge, load_gene_list, zip_files
from datasets import CELL_LINES_TRAINING, CELL_LINES_LEADERBOARD, PRIORITY_GENE_LIST, RESULTS_FOLDER


ESTIMATORS = {'knn': KNeighborsRegressor,  # low scores
              'svm': SVR,
              'lr': LinearRegression,
              'rdg': Ridge,
              'lss': Lasso,  # std = 0
              'eln': ElasticNet,  # std = 0
              'lsscv': LassoCV,  # std = 0
              'elncv': ElasticNetCV,  # std = 0
              'bayr': BayesianRidge,
              'par': PassiveAggressiveRegressor,
              'ard': ARDRegression,  # too slow
              'logr': LogisticRegression,  # too slow
              'gss': GaussianProcess,  # std = 0
              'gssnb': GaussianNB,  # std = 0
              'pls': PLSRegression, #fails
              'tree': DecisionTreeRegressor, # low scores
              'gbr': GradientBoostingRegressor,  #too slow
              'iso': IsotonicRegression  # fails
              }


def average_by_cell_line():
    _, ess_train_data, _ = load_datasets()

    lines_board = load_cell_lines(CELL_LINES_LEADERBOARD)
    lines_train = load_cell_lines(CELL_LINES_TRAINING)

    data = {}

    for line in lines_board.index:
        site = lines_board.at[line, 'Site_primary']
        matches = lines_train.index[lines_train['Site_primary'] == site]
        if matches.size > 0:
            data[line] = ess_train_data.loc[:, matches].mean(1).tolist()
        else:
            data[line] = ess_train_data.mean(1).tolist()

    ess_avg_data = DataFrame(data,
                             index=ess_train_data.index,
                             columns=lines_board.index)

    ess_avg_data.insert(0, 'Description', ess_train_data.index)
    save_gct_data(ess_avg_data, 'avg_per_line.gct')


def pre_process_data(X, Z, method='id', method_args={}):
    X2, Z2 = [], []

    if method == 'id':
        X2, Z2 = X, Z

    if method == 'z-score':
        z_score = X.std(0) / X.mean(0)
        X2 = X[:, z_score > method_args['z_min']]
        Z2 = Z[:, z_score > method_args['z_min']]


    scaler = StandardScaler().fit(X2)
    X2 = scaler.transform(X2)
    Z2 = scaler.transform(Z2)

    return X2, Z2

def pre_process_datasets(train_data, board_data, method='id', method_args={}):

    if method == 'id':
        train_data, board_data

    if method == 'z-score':
        z_score = train_data.std(1).values / train_data.mean(1).values
        train_data = train_data.loc[z_score > method_args['z_min'], :]
        board_data = board_data.loc[z_score > method_args['z_min'], :]


    scaler = StandardScaler().fit(train_data.values.T)
    train_data.values[:,:] = scaler.transform(train_data.values.T).T
    board_data.values[:,:] = scaler.transform(board_data.values.T).T

    return train_data, board_data

def train_and_predict(X, Y, Z, method, method_args):
    W = []
    estimator = ESTIMATORS[method](**method_args)

    for i, y in enumerate(Y):
        estimator.fit(X, y)
        w = estimator.predict(Z)
        W.append(w)
        if (i+1) % (len(Y) / 10) == 0:
            print '.',
    print

    return W


def training_score(Y, Y2):
    return mean([spearmanr(y, y2)[0] for y, y2 in zip(Y, Y2)])


def run_pipeline_sc1(preprocess, method, outputfile, pre_process_args={}, method_args={}, submit=False):
    exp_train_data, ess_train_data, exp_board_data = load_datasets()

    #ess_train_data = ess_train_data.head(100)

    X = exp_train_data.values.T
    Y = ess_train_data.values
    Z = exp_board_data.values.T

    print 'pre-processing data using method', preprocess, str(pre_process_args)
    X2, Z2 = pre_process_data(X, Z, preprocess, pre_process_args)

    t0 = time()
    print 'training and predicting using method', method, str(method_args)
    W = train_and_predict(X2, Y, Z2, method, method_args)
    t1 = time() - t0
    print 'tested', method, str(method_args), 'scored:', 'elapsed', t1, 'secs'

    ess_board_data = DataFrame(W,
                               columns=exp_board_data.columns,
                               index=ess_train_data.index)

    ess_board_data.insert(0, 'Description', ess_train_data.index)
    save_gct_data(ess_board_data, outputfile)

    if submit:
        label = 'daniel_' + outputfile[:-4]
        submit_to_challenge(outputfile, 'sc1', label)


def score_from_training_set(preprocess, method, pre_process_args={}, method_args={}):
    exp_train_data, ess_train_data, _ = load_datasets()

    exp_score_data = exp_train_data.iloc[:,2::3]
    exp_train_data = exp_train_data.iloc[:,::3].join(exp_train_data.iloc[:,1::3])

    ess_train_data = ess_train_data.head(100)

    ess_score_data = ess_train_data.iloc[:,2::3]
    ess_train_data = ess_train_data.iloc[:,::3].join(ess_train_data.iloc[:,1::3])


    X = exp_train_data.values.T
    Y = ess_train_data.values
    Z = exp_score_data.values.T
    W = ess_score_data.values

    print 'pre-processing data using', preprocess, str(pre_process_args)
    X, Z = pre_process_data(X, Z, preprocess, pre_process_args)

    t0 = time()
    print 'training and predicting using', method, str(method_args)
    W2 = train_and_predict(X, Y, Z, method, method_args)
    t1 = time() - t0

    score = training_score(W, W2)

    print 'tested', method, str(method_args), 'scored:', score, 'elapsed', t1, 'secs'


def select_features_per_gene(X, Y, Z, feature_list, selection_method, estimator_method, selection_args, estimator_args):
    W = []
    features = []

    if estimator_method == 'svm':
        estimator_args['kernel'] = 'linear'

    estimator = ESTIMATORS[estimator_method](**estimator_args)

    if selection_method == 'RFE':
        selector = RFE(estimator=estimator, n_features_to_select=10, **selection_args)

    for i, y in enumerate(Y):
        selector = selector.fit(X, y)
        features.append(feature_list[selector.support_])
        w = selector.predict(Z)
        W.append(w)
        if (i+1) % (len(Y) / 10) == 0:
            print '.',
    print


    return W, features


def run_pipeline_sc2(selection_method, estimator_method, outputfile, selection_args={}, estimator_args={}, z_min=0, submit=False):
    exp_train_data, ess_train_data, exp_board_data = load_datasets()
    gene_list = load_gene_list(PRIORITY_GENE_LIST)

    exp_train_data, exp_board_data = pre_process_datasets(exp_train_data, exp_board_data, 'z-score', {'z_min': z_min})

    X = exp_train_data.values.T
    Y = ess_train_data.loc[gene_list, :].values
    Z = exp_board_data.values.T
    feature_list = exp_board_data.index.values


    print 'selection with', selection_method, str(selection_args)
    print 'prediction with', estimator_method, str(estimator_args)

    t0 = time()
    W, features = select_features_per_gene(X, Y, Z, feature_list, selection_method, estimator_method, selection_args, estimator_args)
    t1 = time() - t0
    print 'elasped', t1

    ess_board_data = DataFrame(W,
                               columns=exp_board_data.columns,
                               index=gene_list)
    ess_board_data.index.name = 'Name'

    ess_board_data.insert(0, 'Description', gene_list)
    save_gct_data(ess_board_data, outputfile + '.gct')

    features_data = DataFrame(features, index=gene_list)
    features_data.to_csv(RESULTS_FOLDER + outputfile + '.txt', sep='\t', header=False)

    zip_files(outputfile, [outputfile + '.txt', outputfile + '.gct'])

    if submit:
        label = 'daniel_' + outputfile
        submit_to_challenge(outputfile + '.zip', 'sc2', label)

