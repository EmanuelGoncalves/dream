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
        'sc': 'sc3',                      # sub-challenge
        'filter': 'var',                  # filtering method (None / variance / coefficient of variation)
        'filter_threshold': (0.65, 1),    # filter threshold (gene expression, copy number varation)
        'use_cnv': False,                 # use copy number variation ?
        'use_mut': False,                 # use mutation data ?
        'normalize': False,               # normalize features ?
        'feature_selection': 'kb',        # feature selection (None / KBest (kb) / Recursive Feature Elimination (rfe))
        'n_features': 30,                 # select number of maximum features
        'selection_args': {},             # args to pass to feature selection method
        'estimator': 'woc',               # estimation method (knn, svm, lr, lgr, par, rdg, lss, eln, rdgcv, lsscv, elncv, mtlss, mteln)
        'estimation_args': {},            # args to pass to estimation method
        'submit': False,                  # submit result to challenge ?
        'outputfile': 'out',              # output file name
        'split_train_set': True,          # split training set for score calculation ?
        'max_predictions': None           # limit number of predictions during training (for speed)
    }

    args_list = [update_dict(default_args,
                             {'estimation_args': {'estimators': estimators}})
                 for estimators in [{'nusvm': {}, 'rdgcv': {'normalize': True}, 'par': {'epsilon': 0.01}}]]

    if multi_threaded:
        p = Pool()
        p.map(pipeline, args_list)
    else:
        map(pipeline, args_list)


if __name__ == '__main__':
    run_pipeline()

