__author__ = 'daniel'

from multiprocessing import Pool
from learning import run_pipeline, score_from_training_set


def call_pipeline(args):
    run_pipeline(*args)


def call_score(args):
    score_from_training_set(*args)


def pipeline_single():
    run_pipeline('z-score', 'eln', 'z03_eln_a05_l01.gct', {'z_min': 0.3}, {'alpha': 0.5, 'l1_ratio': 0.1}, 'sc1')


def pipeline_parallel():

    methods = ['lr', 'rdg', 'par']
    z_min = [0.3, 0.4, 0.5]
    args = [('z-score', method, 'z_min_{0}_{1}.gct'.format(z, method), {'z_min': z}, {}, 'sc1')
            for z in z_min for method in methods]
    p = Pool()
    p.map(call_pipeline, args)


def train_set_score():
    score_from_training_set('z-score', 'par', {'z_min': 0.4}, {})
    # score_from_training_set('z-score', 'svm', {'z_min': 0.4}, {})
    # score_from_training_set('z-score', 'ridge', {'z_min': 0.4}, {'alpha': 0.5})


if __name__ == '__main__':
    #average_by_cell_line()
    pipeline_single()
