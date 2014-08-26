__author__ = 'daniel'

from multiprocessing import Pool
from learning import pipeline


def call_pipeline(args):
    pipeline(*args)


def run_pipeline():

    sc = 'sc1'
    preprocess_method = 'cv'
    preprocess_args = {'t': 0.1}
    selection_method = 'KBest'
#    selection_args = {'n_features': 100, 'step': 100}
    selection_args = {'n_features': 2000}
    estimation_method = 'rdgcv'
    estimation_args = {}
    submit = True
    outputfile = 'cv01kb2e3rdgcv'
    multi_threaded = False

    args_list = [(sc, preprocess_method, selection_method, estimation_method, outputfile,
                  preprocess_args, selection_args, estimation_args, submit)]

    if multi_threaded:
        p = Pool()
        p.map(call_pipeline, args_list)
    else:
        map(call_pipeline, args_list)


if __name__ == '__main__':
    #average_by_cell_line()
    run_pipeline()

