__author__ = 'daniel'

from multiprocessing import Pool
from learning import run_pipeline

def call_pipeline(args):
    run_pipeline(*args)

def main():
    #average_by_cell_line()

    methods = ['bayr', 'par']
    z_min = [0.1, 0.2]
    args = [('z-score', method, 'z_min_{0}_{1}.gct'.format(z, method), {'z_min': z}, {})
            for method in methods for z in z_min]
    p = Pool()
    p.map(call_pipeline, args)


if __name__ == '__main__':
    main()
