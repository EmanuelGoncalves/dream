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
        'sc': 'sc1',
        'filter': 'cv',
        'filter_threshold': (0.1, 0.5),
        'use_cnv': False,
        'normalize': True,
        'feature_selection': None,
        'n_features': 5000,
        'selection_args': {},
        'estimator': 'rdg',
        'estimation_args': {},
        'submit': True,
        'outputfile': 'out',
        'split_train_set': False
    }

    args_list = [update_dict(default_args,
                             {'outputfile': 'cn{0}_norm{1}_cv01-05_{2}'.format(use_cnv, normalize, method),
                              'use_cnv': bool(use_cnv),
                              'normalize': bool(normalize)})
                 for normalize in [0, 1]
                 for use_cnv in [0, 1]
                 for method in ['par', 'rdg']]

    if multi_threaded:
        p = Pool()
        p.map(pipeline, args_list)
    else:
        map(pipeline, args_list)


if __name__ == '__main__':
    run_pipeline()

