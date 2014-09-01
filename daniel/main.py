__author__ = 'daniel'

from multiprocessing import Pool
from learning import pipeline


def update_dict(x, y):
    z = x.copy()
    z.update(y)
    return z

def run_pipeline():

    multi_threaded = False

    default_args = {
        'sc': 'sc3',
        'filter': 'var',
        'filter_threshold': (0, 0),
        'use_cnv': False,
        'normalize': False,
        'feature_selection': 'kb',
        'n_features': 100,
        'selection_args': {},
        'estimator': 'rdgcv',
        'estimation_args': {},
        'submit': True,
        'outputfile': 'out',
        'split_train_set': False,
        'max_predictions': 1000
    }

    args_list = [update_dict(default_args,
                             {'outputfile': 'best_of_emanuel'})]

    if multi_threaded:
        p = Pool()
        p.map(pipeline, args_list)
    else:
        map(pipeline, args_list)


if __name__ == '__main__':
    run_pipeline()

