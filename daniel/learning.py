__author__ = 'daniel'

from pandas import DataFrame
from numpy import mean, std
from numpy.random import randn
from scipy.stats import spearmanr
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet, LassoCV, \
    ElasticNetCV, PassiveAggressiveRegressor, MultiTaskLasso, RidgeCV
from sklearn.feature_selection import RFE, SelectKBest, f_regression, VarianceThreshold
from time import time

from datasets import load_datasets, load_cell_lines, save_gct_data, submit_to_challenge, load_gene_list, zip_files
from datasets import CELL_LINES_TRAINING, CELL_LINES_LEADERBOARD, PRIORITY_GENE_LIST, RESULTS_FOLDER


ESTIMATORS = {'knn': KNeighborsRegressor,
              'svm': SVR,
              'lr': LinearRegression,
              'rdg': Ridge,
              'rdgcv': RidgeCV,
              'lss': Lasso,
              'lsscv': LassoCV,
              'mlss': MultiTaskLasso,
              'eln': ElasticNet,
              'elncv': ElasticNetCV,
              'par': PassiveAggressiveRegressor,
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

    ess_avg_data = DataFrame(data, index=ess_train_data.index, columns=lines_board.index)
    ess_avg_data.insert(0, 'Description', ess_train_data.index)
    save_gct_data(ess_avg_data, 'avg_per_line.gct')


def pre_process_datasets(train_data, board_data, filter_method=None, method_args={}):

    if filter_method is None:
        train_data, board_data

    if filter_method == 'cv':
        cv = train_data.std(1).values / train_data.mean(1).values
        train_data = train_data.loc[cv > method_args['t'], :]
        board_data = board_data.loc[cv > method_args['t'], :]

    if filter_method == 'var':
        selector = VarianceThreshold(method_args['t'])
        selector.fit(train_data.values.T)
        train_data = train_data.loc[selector.get_support(), :]
        board_data = board_data.loc[selector.get_support(), :]

    scaler = StandardScaler().fit(train_data.values.T)
    train_data.values[:,:] = scaler.transform(train_data.values.T).T
    board_data.values[:,:] = scaler.transform(board_data.values.T).T

    return train_data, board_data


def select_train_predict(X, Y, Z, feature_list, selection_method, estimator_method, selection_args, estimator_args):
    W = []
    features = []

    if estimator_method == 'svm' and selection_method == 'RFE':
        estimator_args['kernel'] = 'linear'

    estimator = ESTIMATORS[estimator_method](**estimator_args)

    n_features = min(len(feature_list), selection_args['n_features'])

    if selection_method is None:
        for i, y in enumerate(Y):
            estimator.fit(X, y)
            w = estimator.predict(Z)
            W.append(w)
            if (i+1) % (len(Y) / 10) == 0:
                print '.',

    if selection_method == 'RFE':
        selector = RFE(estimator=estimator, n_features_to_select=n_features)

        for i, y in enumerate(Y):
            selector = selector.fit(X, y)
            features.append(feature_list[selector.support_])
            w = selector.predict(Z)
            W.append(w)
            if (i+1) % (len(Y) / 10) == 0:
                print '.',

    if selection_method == 'KBest':
        selector = SelectKBest(f_regression, k=n_features)
        for i, y in enumerate(Y):
            X2 = selector.fit_transform(X, y)
            Z2 = selector.transform(Z)
            features.append(feature_list[selector.get_support()])
            estimator.fit(X2, y)
            w = estimator.predict(Z2)
            W.append(w)
            if (i+1) % (len(Y) / 10) == 0:
                print '.',

    print

    bad = 0
    for i, w in enumerate(W):
        if std(w) < 1e-12:
            W[i] = randn(*w.shape)
            bad += 1

    if bad:
        print '* Warning:', bad, 'bad predictions out of', len(W)

    return W, features


def sc3_multitask(X, Y, Z, feature_list, selection_method, estimator_method, selection_args, estimator_args):

    W = []
    features = []

    if estimator_method == 'svm' and selection_method == 'RFE':
        estimator_args['kernel'] = 'linear'

    n_features = min(len(feature_list), selection_args['n_features'])

    estimator = ESTIMATORS[estimator_method](**estimator_args)

    if selection_method == 'RFE':
        del selection_args['n_features']
        selector = RFE(estimator=estimator, n_features_to_select=n_features, **selection_args)
        selector = selector.fit(X, Y.T)
        features = feature_list[selector.support_]
        W = selector.predict(Z)

    if selection_method == 'KBest':
        print 'Cannot use KBest with multi task methods'


    return W.T, features



def sc3_top100(X, Y, Z, feature_list, selection_method, estimation_method, selection_args, estimation_args):

    _, features = select_train_predict(X, Y, Z, feature_list, selection_method, estimation_method, selection_args, estimation_args)

    best_features = [feat for row in features for feat in row]
    feature_count = [(feat, best_features.count(feat)) for feat in set(best_features)]
    feature_count.sort(key=lambda (a, b): b, reverse=True)
    top100 = [feat for feat, _ in feature_count[:100]]
    feature_list = feature_list.tolist()
    indices = [feature_list.index(item) for item in top100]
    X2 = X[:,indices]
    Z2 = Z[:,indices]

    estimator = ESTIMATORS[estimation_method](**estimation_args)
    W = []

    for i, y in enumerate(Y):
        estimator.fit(X2, y)
        w = estimator.predict(Z2)
        W.append(w)
        if (i+1) % (len(Y) / 10) == 0:
            print '.',
    print

    return W, top100


def pipeline(sc, preprocess_method, selection_method, estimation_method, outputfile, preprocess_args={}, selection_args={}, estimation_args={}, submit=False):
    exp_train_data, ess_train_data, exp_board_data = load_datasets()
    gene_list = load_gene_list(PRIORITY_GENE_LIST)

    print 'pre-processing with:', preprocess_method, str(preprocess_args)

    exp_train_data, exp_board_data = pre_process_datasets(exp_train_data, exp_board_data, preprocess_method, preprocess_args)

    X = exp_train_data.values.T
    Y = ess_train_data.values if sc == 'sc1' else ess_train_data.loc[gene_list, :].values
    Z = exp_board_data.values.T
    feature_list = exp_board_data.index.values
    print 'features after filtering:', len(feature_list)

    print 'predicting with', selection_method, estimation_method,

    t0 = time()
    if sc == 'sc1':
        W, features = select_train_predict(X, Y, Z, feature_list, selection_method, estimation_method, selection_args, estimation_args)

    if sc == 'sc2':
        selection_args['n_features'] == 10
        W, features = select_train_predict(X, Y, Z, feature_list, selection_method, estimation_method, selection_args, estimation_args)

    if sc == 'sc3':
        if estimation_method == 'mlss':
            selection_args['n_features'] == 100
            W, features = sc3_multitask(X, Y, Z, feature_list, selection_method, estimation_method, selection_args, estimation_args)
        else:
            W, features = sc3_top100(X, Y, Z, feature_list, selection_method, estimation_method, selection_args, estimation_args)

    t1 = time() - t0
    print 'tested', selection_method, estimation_method, 'elapsed', t1, 'secs'

    index = ess_train_data.index if sc == 'sc1' else gene_list
    ess_board_data = DataFrame(W, columns=exp_board_data.columns, index=index)
    ess_board_data.insert(0, 'Description', index)

    if sc == 'sc1':
        save_gct_data(ess_board_data, outputfile + '.gct')

    if sc == 'sc2':
        features_data = DataFrame(features, index=gene_list)
        features_data.to_csv(RESULTS_FOLDER + outputfile + '.txt', sep='\t', header=False)
        zip_files(outputfile + '.zip', [outputfile + '.txt', outputfile + '.gct'])

    if sc == 'sc3':
        features_data = DataFrame([features])
        features_data.to_csv(RESULTS_FOLDER + outputfile + '.txt', sep='\t', header=False, index=False)
        zip_files(outputfile + '.zip', [outputfile + '.txt', outputfile + '.gct'])

    if submit:
        label = 'ph2_daniel_' + outputfile[:-4]
        submit_to_challenge(outputfile, sc, label)


def training_score(Y, Y2):
    return mean([spearmanr(y, y2)[0] for y, y2 in zip(Y, Y2)])



def score_from_training_set(preprocess, method, pre_process_args={}, method_args={}, limit_to=None):
    exp_train_data, ess_train_data, _ = load_datasets()


    exp_score_data = exp_train_data.iloc[:,2::3]
    exp_train_data = exp_train_data.iloc[:,::3].join(exp_train_data.iloc[:,1::3])

    print 'pre-processing data using method', preprocess, str(pre_process_args)
    exp_train_data, exp_score_data = pre_process_datasets(exp_train_data, exp_score_data, preprocess, pre_process_args)


    if limit_to:
        ess_train_data = ess_train_data.head(limit_to)

    ess_score_data = ess_train_data.iloc[:,2::3]
    ess_train_data = ess_train_data.iloc[:,::3].join(ess_train_data.iloc[:,1::3])


    X = exp_train_data.values.T
    Y = ess_train_data.values
    Z = exp_score_data.values.T
    W = ess_score_data.values

    t0 = time()
    print 'training and predicting using', method, str(method_args)
    W2 = select_train_predict(X, Y, Z, method, method_args)
    t1 = time() - t0

    score = training_score(W, W2)

    print 'tested', method, str(method_args), 'scored:', score, 'elapsed', t1, 'secs'
