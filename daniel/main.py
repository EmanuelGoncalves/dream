__author__ = 'daniel'

from multiprocessing import Pool
from learning import run_pipeline

def call_pipeline(args):
    run_pipeline(*args)

def main():
    #average_by_cell_line()

    methods = ['ridge', 'lasso', 'elnet', 'linreg']
    args = [('z-score', method, 'z_min_01_ridge.gct', {'z_min': 0.08}, {})  for method in methods]
    p = Pool()
    p.map(call_pipeline, args)


if __name__ == '__main__':
    main()
