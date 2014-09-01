__author__ = 'daniel'

from multiprocessing import Pool
from learning import pipeline


def update_dict(x, y):
    z = x.copy()
    z.update(y)
    return z

def run_pipeline():

    multi_threaded = True

    default_args = {
        'sc': 'sc2',
        'filter': 'var',
        'filter_threshold': (0.65, 1),
        'use_cnv': False,
        'normalize': False,
        'feature_selection': 'kb',
        'n_features': 1000,
        'selection_args': {},
        'estimator': 'rdgcv',
        'estimation_args': {},
        'submit': True,
        'outputfile': 'out',
        'split_train_set': False,
        'max_predictions': 100
    }

    args_list = [update_dict(default_args,
                             {'outputfile': 'n{0}_cnv{1}_{2}'.format(norm, cnv, method),
                             'normalize': bool(norm),
                             'use_cnv': bool(cnv),
                             'estimator': method})
                 for method in ['rdgcv', 'svm']
                 for norm in [0, 1]
                 for cnv in [0, 1]]

    if multi_threaded:
        p = Pool()
        p.map(pipeline, args_list)
    else:
        map(pipeline, args_list)


if __name__ == '__main__':
    run_pipeline()

