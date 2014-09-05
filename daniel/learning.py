__author__ = 'daniel'

from pandas import DataFrame
from numpy import mean, std, exp
from numpy.random import randn
from scipy.stats import spearmanr
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression, LogisticRegression, Ridge, Lasso, ElasticNet, PassiveAggressiveRegressor, \
    RidgeCV, LassoCV, ElasticNetCV, MultiTaskLasso, MultiTaskElasticNet
from sklearn.feature_selection import RFE, SelectKBest, f_regression, VarianceThreshold
from time import time

from datasets import load_datasets, save_gct_data, submit_to_challenge, load_gene_list, zip_files
from datasets import RESULTS_FOLDER


ESTIMATORS = {'knn': KNeighborsRegressor,
              'svm': SVR,
              'lr': LinearRegression,
              'lgr': LogisticRegression,
              'par': PassiveAggressiveRegressor,
              'rdg': Ridge,
              'lss': Lasso,
              'eln': ElasticNet,
              'rdgcv': RidgeCV,
              'lsscv': LassoCV,
              'elncv': ElasticNetCV,
              'mtlss': MultiTaskLasso,
              'mteln': MultiTaskElasticNet,
              }


def pre_process_datasets(datasets, filter_method=None, threshold=(0, 0), normalize=True, use_cnv=False, use_mut=False):

    exp_train_data = datasets['exp_train_data']
    exp_board_data = datasets['exp_board_data']

    if use_cnv:
        cnv_train_data = datasets['cnv_train_data']
        cnv_board_data = datasets['cnv_board_data']

    if filter_method == 'cv':
        exp_cv = exp_train_data.std(1).values / exp_train_data.mean(1).values
        exp_train_data = exp_train_data.loc[exp_cv > threshold[0], :]
        exp_board_data = exp_board_data.loc[exp_cv > threshold[0], :]
        if use_cnv:
            cnv_train_data = cnv_train_data.apply(exp)
            cnv_cv = cnv_train_data.std(1).values / cnv_train_data.mean(1).values
            cnv_train_data = cnv_train_data.loc[cnv_cv > threshold[1], :]
            cnv_board_data = cnv_board_data.loc[cnv_cv > threshold[1], :]

    if filter_method == 'var':
        selector = VarianceThreshold(threshold[0])
        selector.fit(exp_train_data.values.T)
        exp_train_data = exp_train_data.loc[selector.get_support(), :]
        exp_board_data = exp_board_data.loc[selector.get_support(), :]
        if use_cnv:
            selector = VarianceThreshold(threshold[1])
            selector.fit(cnv_train_data.values.T)
            cnv_train_data = cnv_train_data.loc[selector.get_support(), :]
            cnv_board_data = cnv_board_data.loc[selector.get_support(), :]

    if use_cnv:
        feat_train_data = exp_train_data.append(cnv_train_data)
        feat_board_data = exp_board_data.append(cnv_board_data)
        print 'features after filtering', exp_train_data.shape[0], '+', cnv_train_data.shape[0], '=', feat_train_data.shape[0]
    else:
        feat_train_data = exp_train_data
        feat_board_data = exp_board_data
        print 'features after filtering', exp_train_data.shape[0]

    if normalize:
        scaler = StandardScaler().fit(feat_train_data.values.T)
        feat_train_data.values[:,:] = scaler.transform(feat_train_data.values.T).T
        feat_board_data.values[:,:] = scaler.transform(feat_board_data.values.T).T

    if use_mut:
        feat_train_data = feat_train_data.append(datasets['mut_train_data'])
        feat_board_data = feat_board_data.append(datasets['mut_board_data'])

    datasets['feat_train_data'] = feat_train_data
    datasets['feat_board_data'] = feat_board_data



def select_train_predict(X, Y, Z, feature_list, selection_method, estimator_method, n_features, selection_args, estimator_args):
    W = []
    features = []
    n_features = min(n_features, len(feature_list))

    if estimator_method == 'svm' and selection_method == 'rfe':
        estimator_args['kernel'] = 'linear'

    estimator = ESTIMATORS[estimator_method](**estimator_args)

    if selection_method is None:
        for i, y in enumerate(Y):
            estimator.fit(X, y)
            w = estimator.predict(Z)
            W.append(w)
            if (i+1) % (len(Y) / 10) == 0:
                print '.',

    if selection_method == 'rfe':
        selector = RFE(estimator=estimator, n_features_to_select=n_features, **selection_args)

        for i, y in enumerate(Y):
            selector = selector.fit(X, y)
            features.append(feature_list[selector.support_])
            w = selector.predict(Z)
            W.append(w)
            if (i+1) % (len(Y) / 10) == 0:
                print '.',

    if selection_method == 'kb':
        selector = SelectKBest(f_regression, k=n_features, **selection_args)
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

    return W, features


def sc3_multitask(X, Y, Z, feature_list, selection_method, estimator_method, selection_args, estimator_args):

    W = []
    features = []

    if estimator_method == 'svm' and selection_method == 'RFE':
        estimator_args['kernel'] = 'linear'

    n_features = min(len(feature_list), selection_args['n_features'])

    estimator = ESTIMATORS[estimator_method](**estimator_args)

    if selection_method == 'rfe':
        del selection_args['n_features']
        selector = RFE(estimator=estimator, n_features_to_select=n_features, **selection_args)
        selector = selector.fit(X, Y.T)
        features = feature_list[selector.support_]
        W = selector.predict(Z)

    if selection_method == 'kb':
        print 'Cannot use KBest with multi task methods'

    return W.T, features


def sc3_top100(X, Y, Z, feature_list, selection_method, estimation_method, n_features, selection_args, estimation_args):

    _, features = select_train_predict(X, Y, Z, feature_list, selection_method, estimation_method, n_features, selection_args, estimation_args)

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


def clean_bad_predictions(W):

    bad = 0
    for i, w in enumerate(W):
        if std(w) < 1e-12:
            W[i] = randn(*w.shape)
            bad += 1

    if bad:
        print '* Warning:', bad, 'bad predictions out of', len(W)


def training_score(X, Y):
    return mean([spearmanr(x, y)[0] for x, y in zip(X, Y)])


def split_datasets(datasets, use_cnv, use_mut):

    idx_train = [0,   1,   2,   3,   8,   9,  10,  12,  14,  16,  21,  23,  24,
                 25,  27,  28,  29,  30,  32,  34,  36,  37,  38,  39,  40,  42,
                 43,  44,  45,  47,  48,  49,  52,  53,  56,  57,  58,  59,  60,
                 61,  62,  63,  64,  66,  67,  68,  69,  70,  73,  75,  80,  82,
                 83,  84,  88,  90,  91,  92,  94,  95,  96,  98,  99, 100, 102, 103]

    idx_board = [4,   5,   6,   7,  11,  13,  15,  17,  18,  19,  20,  22,  26,
                 31,  33,  35,  41,  46,  50,  51,  54,  55,  65,  71,  72,  74,
                 76,  77,  78,  79,  81,  85,  86,  87,  89,  93,  97, 101, 104]

    exp_train_data = datasets['exp_train_data'][idx_train]
    datasets['exp_board_data'] = datasets['exp_train_data'][idx_board]
    datasets['exp_train_data'] = exp_train_data

    ess_train_data = datasets['ess_train_data'][idx_train]
    datasets['ess_score_data'] = datasets['ess_train_data'][idx_board]
    datasets['ess_train_data'] = ess_train_data

    if use_cnv:
        cnv_train_data = datasets['cnv_train_data'][idx_train]
        datasets['cnv_board_data'] = datasets['cnv_train_data'][idx_board]
        datasets['cnv_train_data'] = cnv_train_data

    if use_mut:
        mut_train_data = datasets['mut_train_data'][idx_train]
        datasets['mut_board_data'] = datasets['mut_train_data'][idx_board]
        datasets['mut_train_data'] = mut_train_data


def pipeline(args):
    phase = args['phase']
    sc = args['sc']
    filter_method = args['filter']
    filter_threshold = args['filter_threshold']
    normalize = args['normalize']
    feature_selection = args['feature_selection']
    n_features = args['n_features'] if sc != 'sc2' else 10
    selection_args = args['selection_args']
    estimator = args['estimator']
    estimation_args = args['estimation_args']
    submit = args['submit']
    outputfile = args['outputfile']
    use_cnv = args['use_cnv']
    use_mut = args['use_mut'] if phase == 'phase3' else False
    split_train_set = args['split_train_set']
    max_predictions = args['max_predictions']

    datasets = load_datasets(phase=phase, get_cnv=use_cnv, get_mut=use_mut)
    gene_list = datasets['gene_list']

    if split_train_set:
        split_datasets(datasets, use_cnv, use_mut)

    print 'pre-processing with:', filter_method, 'at', filter_threshold, 'normalize:', normalize, 'use_cnv', use_cnv, 'use_mut', use_mut

    pre_process_datasets(datasets, filter_method, filter_threshold, normalize, use_cnv)

    feat_train_data = datasets['feat_train_data']
    feat_board_data = datasets['feat_board_data']
    ess_train_data = datasets['ess_train_data']
    X = feat_train_data.values.T
    Y = ess_train_data.values if sc == 'sc1' else ess_train_data.loc[gene_list, :].values
    Z = feat_board_data.values.T
    feature_list = feat_board_data.index.values

    if max_predictions and split_train_set:
        Y = Y[::len(Y)/max_predictions][:max_predictions]

    print 'predicting with', feature_selection, '(', n_features, ')',  estimator,

    t0 = time()
    if sc == 'sc1':
        W, features = select_train_predict(X, Y, Z, feature_list, feature_selection, estimator, n_features, selection_args, estimation_args)

    if sc == 'sc2':
        W, features = select_train_predict(X, Y, Z, feature_list, feature_selection, estimator, n_features, selection_args, estimation_args)

    if sc == 'sc3':
        if estimator.startswith('mt'):
            W, features = sc3_multitask(X, Y, Z, feature_list, feature_selection, estimator, selection_args, estimation_args)
        else:
            W, features = sc3_top100(X, Y, Z, feature_list, feature_selection, estimator, n_features, selection_args, estimation_args)

    t1 = time() - t0
    clean_bad_predictions(W)

    print 'tested', feature_selection,  estimator, 'elapsed', t1, 'secs'

    if split_train_set:
        W0 = datasets['ess_score_data'].values if sc == 'sc1' else datasets['ess_score_data'].loc[gene_list, :].values
        if max_predictions:
            W0 = W0[::len(W0)/max_predictions][:max_predictions]
        score = training_score(W0, W)
        print 'scored:', score

    else:
        index = ess_train_data.index if sc == 'sc1' else gene_list
        ess_board_data = DataFrame(W, columns=feat_board_data.columns, index=index)
        ess_board_data.insert(0, 'Description', index)
        save_gct_data(ess_board_data, outputfile + '.gct')

        if sc == 'sc2':
            features_data = DataFrame(features, index=gene_list)
            features_data.to_csv(RESULTS_FOLDER + outputfile + '.txt', sep='\t', header=False)
            zip_files(outputfile, [outputfile + '.txt', outputfile + '.gct'])

        if sc == 'sc3':
            features_data = DataFrame([features])
            features_data.to_csv(RESULTS_FOLDER + outputfile + '.txt', sep='\t', header=False, index=False)
            zip_files(outputfile, [outputfile + '.txt', outputfile + '.gct'])

        if submit:
            label = 'daniel_' + outputfile
            submit_to_challenge(outputfile, sc, label)

