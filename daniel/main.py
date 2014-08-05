__author__ = 'daniel'

from learning import run_pipeline

def main():
    #average_by_cell_line()
    #run_pipeline('z-score', 'knn', 'z_min_02_knn_5.gct', {'z_min': 0.2}, {'n_neighbors': 5})
    #run_pipeline('z-score', 'svm', 'z_min_005_svm.gct', {'z_min': 0.05}, {})
    #run_pipeline('z-score', 'svm', 'z_min_008_svm.gct', {'z_min': 0.08}, {})
    #run_pipeline('z-score', 'svm', 'z_min_015_svm.gct', {'z_min': 0.15}, {})
    #run_pipeline('z-score', 'linreg', 'z_min_02_linreg.gct', {'z_min': 0.2}, {})
    run_pipeline('z-score', 'linreg', 'z_min_01_linreg.gct', {'z_min': 0.1}, {})


if __name__ == '__main__':
    main()
    