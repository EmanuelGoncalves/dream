__author__ = 'daniel'

from pandas import DataFrame
from numpy import mean, std, exp
from numpy.random import randn
from scipy.stats import spearmanr
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet, PassiveAggressiveRegressor, \
    RidgeCV, LassoCV, ElasticNetCV, MultiTaskLasso, MultiTaskElasticNet
from sklearn.feature_selection import RFE, SelectKBest, f_regression, VarianceThreshold
from time import time

from datasets import load_datasets, save_gct_data, submit_to_challenge, load_gene_list, zip_files
from datasets import RESULTS_FOLDER


ESTIMATORS = {'knn': KNeighborsRegressor,
              'svm': SVR,
              'lr': LinearRegression,
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


def pre_process_datasets(datasets, filter_method=None, threshold=(0, 0), normalize=True, use_cnv=False):

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

    feat_train_data = exp_train_data.append(cnv_train_data) if use_cnv else exp_train_data
    feat_board_data = exp_board_data.append(cnv_board_data) if use_cnv else exp_board_data

    if normalize:
        scaler = StandardScaler().fit(feat_train_data.values.T)
        feat_train_data.values[:,:] = scaler.transform(feat_train_data.values.T).T
        feat_board_data.values[:,:] = scaler.transform(feat_board_data.values.T).T

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


def split_datasets(datasets, use_cnv):
    exp_train_data = datasets['exp_train_data'].iloc[:,::3].join(datasets['exp_train_data'].iloc[:,1::3])
    datasets['exp_board_data'] = datasets['exp_train_data'].iloc[:,2::3]
    datasets['exp_train_data'] = exp_train_data

    ess_train_data = datasets['ess_train_data'].iloc[:,::3].join(datasets['ess_train_data'].iloc[:,1::3])
    datasets['ess_score_data'] = datasets['ess_train_data'].iloc[:,2::3]
    datasets['ess_train_data'] = ess_train_data

    if use_cnv:
        cnv_train_data = datasets['cnv_train_data'].iloc[:,::3].join(datasets['cnv_train_data'].iloc[:,1::3])
        datasets['cnv_board_data'] = datasets['cnv_train_data'].iloc[:,2::3]
        datasets['cnv_train_data'] = cnv_train_data


def pipeline(args):
    sc = args['sc']
    filter_method = args['filter']
    filter_threshold = args['filter_threshold']
    normalize = args['normalize']
    feature_selection = args['feature_selection']
    n_features = args['n_features']
    selection_args = args['selection_args']
    estimator = args['estimator']
    estimation_args = args['estimation_args']
    submit = args['submit']
    outputfile = args['outputfile']
    use_cnv = args['use_cnv']
    split_train_set = args['split_train_set']

    datasets = load_datasets(get_cnv=use_cnv)
    gene_list = datasets['gene_list']

    if split_train_set:
        split_datasets(datasets, use_cnv)

    print 'pre-processing with:', filter_method, 'at', filter_threshold, 'normalize:', normalize

    pre_process_datasets(datasets, filter_method, filter_threshold, normalize, use_cnv)

    feat_train_data = datasets['feat_train_data']
    feat_board_data = datasets['feat_board_data']
    ess_train_data = datasets['ess_train_data']
    X = feat_train_data.values.T
    Y = ess_train_data.values if sc == 'sc1' else ess_train_data.loc[gene_list, :].values
    Z = feat_board_data.values.T
    feature_list = feat_board_data.index.values

    print 'features after filtering:', len(feature_list)

    print 'predicting with', feature_selection,  estimator,

    t0 = time()
    if sc == 'sc1':
        W, features = select_train_predict(X, Y, Z, feature_list, feature_selection, estimator, n_features, selection_args, estimation_args)

    if sc == 'sc2':
        n_features = 10
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
        W0 = datasets['ess_score_data'].values
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

