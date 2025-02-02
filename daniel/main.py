__author__ = 'daniel'

from multiprocessing import Pool
from learning import pipeline


def update_dict(x, y):
    z = x.copy()
    z.update(y)
    return z


def run_pipeline():

    multi_threaded = False                # run multiple processes per configuration ?

    default_args = {
        'phase': 'phase3',                # challenge phase
        'sc': 'sc1',                      # sub-challenge
        'filter': 'var',                  # filtering method (None / variance / coefficient of variation)
        'filter_threshold': (0.65, 0.3),  # filter threshold (gene expression, copy number varation)
        'use_cnv': True,                  # use copy number variation ?
        'use_mut': False,                 # use mutation data ?
        'normalize': False,               # normalize features ?
        'feature_selection': 'kbest',     # feature selection (None / KBest (kb) / Recursive Feature Elimination (rfe))
        'n_features': 3500,               # select number of maximum features
        'selection_args': {},             # args to pass to feature selection method
        'estimator': 'woc',               # estimation method (knn, svm, lr, lgr, par, rdg, lss, eln, rdgcv, lsscv, elncv, mtlss, mteln)
        'estimation_args': {},            # args to pass to estimation method
        'submit': True,                   # submit result to challenge ?
        'outputfile': 'sc1_daniel_phase3',# output file name
        'split_train_set': False,         # split training set for score calculation ?
        'max_predictions': None,          # limit number of predictions during training (for speed)
    }

    args_list = [update_dict(default_args,
                             {'estimation_args': {'estimators': estimators}})
                 for estimators in [{'par': {'epsilon': 0.01}, 'rdgcv': {'normalize': True}, 'nusvm':{}}]]

    if multi_threaded:
        p = Pool()
        p.map(pipeline, args_list)
    else:
        map(pipeline, args_list)


if __name__ == '__main__':
    run_pipeline()

