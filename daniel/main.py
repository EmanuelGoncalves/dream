__author__ = 'daniel'

from multiprocessing import Pool
from learning import run_pipeline_sc1, score_from_training_set, run_pipeline_sc2


def run_sc1(args):
    run_pipeline_sc1(*args)

def run_sc2(args):
    run_pipeline_sc2(*args)


def call_score(args):
    score_from_training_set(*args)


def pipeline_sc1_single():
    run_pipeline_sc1('z-score', 'eln', 'z03_eln_a05_l01.gct', {'z_min': 0.3}, {'alpha': 0.5, 'l1_ratio': 0.1}, True)


def pipeline_sc1_parallel():

    methods = ['lr', 'rdg', 'par']
    z_min = [0.3, 0.4, 0.5]
    args = [('z-score', method, 'z_min_{0}_{1}.gct'.format(z, method), {'z_min': z}, {}, True)
            for z in z_min for method in methods]
    p = Pool()
    p.map(run_sc1, args)


def train_set_score():
    score_from_training_set('z-score', 'elncv', {'z_min': 0.3}, {'l1_ratio': 0.1, 'alpha': 0.1})


def train_set_score_parallel():
    methods = ['par']
    z_min = [0.3]
    C = [0.1, 0.5, 1.0, 1.5, 2.0]
    args = [('z-score', method, {'z_min': z}, {'C': c})
            for z in z_min for method in methods for c in C]
    p = Pool()
    p.map(call_score, args)

def pipeline_sc2_single():
#    run_pipeline_sc2('RFE', 'svm', 'rfe_svm', {'step': 100}, {'kernel': 'linear'}, True)
    run_pipeline_sc2('RFE', 'lr', 'rfe_lr', {'step': 50}, {}, z_min=0.4, submit=True)
    run_pipeline_sc2('RFE', 'par', 'rfe_par', {'step': 50}, {}, z_min=0.4, submit=True)

def pipeline_sc2_parallel():
    methods = ['svm', 'lr', 'rdg', 'par']
    z_min = [0.2, 0.1]
    steps = [100, 50, 10]
    args = [('RFE', method, 'rfe_{0}_z{1}_st{2}'.format(method, z, step), {'step': step}, {}, z, True)
            for z in z_min for step in steps for method in methods]
    p = Pool()
    p.map(run_sc2, args)


if __name__ == '__main__':
    #average_by_cell_line()
    #pipeline_single()
    #pipeline_parallel()
    #train_set_score()
    #train_set_score_parallel()
    #pipeline_sc2_single()
    pipeline_sc2_parallel()
