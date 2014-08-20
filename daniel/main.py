__author__ = 'daniel'

from multiprocessing import Pool
from learning import run_pipeline_sc1, score_from_training_set, run_pipeline_sc2, run_pipeline_sc3


def run_sc1(args):
    run_pipeline_sc1(*args)

def run_sc2(args):
    run_pipeline_sc2(*args)


def call_score(args):
    score_from_training_set(*args)


def pipeline_sc1_single():
    run_pipeline_sc1('var', 'par', 'v_eln_a05_l01.gct', {'t': 0.25}, {}, False)


def pipeline_sc1_parallel():

    methods = ['par']
    T = [0.1, 0.2, 0.3, 0.4, 0.5]
    args = [('var', method, 'v{0}_{1}.gct'.format(t, method), {'t': t}, {}, True)
            for t in T for method in methods]
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
#    run_pipeline_sc2('RFE', 'lr', 'rfe_lr', {'step': 50}, {}, z_min=0.4, submit=True)
    run_pipeline_sc2('KBest', 'rdgcv', 'sc2_t0rdgcv', {}, {}, t=0, submit=True)

def pipeline_sc2_parallel():
    methods = ['lr', 'rdg', 'par']
    args = [('KBest', method, 'kb'+method, {}, {}, 0, True)
            for method in methods]
    p = Pool()
    p.map(run_sc2, args)


def pipeline_sc3_single():
#    run_pipeline_sc3('RFE', 'mlss', 'sc3_rfe_mlss', {'step': 100}, {'alpha': 0.1}, z_min=0.3, submit=True)
    run_pipeline_sc3('KBest', 'rdgcv', 'sc3_kb_rdgcv', {}, {}, z_min=0, submit=False)


if __name__ == '__main__':
    #average_by_cell_line()
    #pipeline_sc1_single()
    #pipeline_sc1_parallel()
    #train_set_score()
    #train_set_score_parallel()
    #pipeline_sc2_single()
    #pipeline_sc2_parallel()
    pipeline_sc3_single()
