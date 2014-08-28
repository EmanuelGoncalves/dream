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
        'filter_threshold': (0.1, 0),
        'use_cnv': False,
        'normalize': True,
        'feature_selection': None,
        'n_features': 5000,
        'selection_args': {},
        'estimator': 'lr',
        'estimation_args': {},
        'submit': True,
        'outputfile': 'out'
    }

    args_list = [update_dict(default_args,
                             {'outputfile': 'norm_{0}'.format(normalize),
                              'normalize': normalize})
                 for normalize in [True, False]]

    if multi_threaded:
        p = Pool()
        p.map(pipeline, args_list)
    else:
        map(pipeline, args_list)


if __name__ == '__main__':
    run_pipeline()

