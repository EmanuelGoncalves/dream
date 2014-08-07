__author__ = 'daniel'

from pandas import read_csv

SC1_DATA = '../data/'

RESULTS_FOLDER = '../submissions/'

CNV_TRAINING_DATA = SC1_DATA + 'CCLE_copynumber_training.gct'
EXP_TRAINING_DATA = SC1_DATA + 'CCLE_expression_training.gct'
ESS_TRAINING_DATA = SC1_DATA + 'Achilles_v2.9_training.gct'
CELL_LINES_TRAINING = SC1_DATA + 'Cell_line_annotation_training.txt'

CNV_LEADERBOARD_DATA = SC1_DATA + 'CCLE_copynumber_leaderboard.gct'
EXP_LEADERBOARD_DATA = SC1_DATA + 'CCLE_expression_leaderboard.gct'
CELL_LINES_LEADERBOARD = SC1_DATA + 'Cell_line_annotation_leaderboard.txt'


def load_gct_data(filename):
    data = read_csv(filename, sep='\t', header=2, index_col=0)
    return data


def save_gct_data(dataframe, filename, folder=RESULTS_FOLDER):
    dataframe.to_csv(folder + 'temp.csv', sep='\t')
    with open(folder + 'temp.csv', 'r') as file1:
        with open(folder + filename, 'w') as file2:
            file2.write('#1.2\n')
            file2.write('{0}\t{1}\n'.format(dataframe.shape[0], dataframe.shape[1] - 1))
            for line in file1:
                file2.write(line)


def load_cell_lines(filename):
    data = read_csv(filename, sep='\t', header=0, index_col=0)
    return data


def load_datasets():
    exp_train_data = load_gct_data(EXP_TRAINING_DATA)
    ess_train_data = load_gct_data(ESS_TRAINING_DATA)
    exp_board_data = load_gct_data(EXP_LEADERBOARD_DATA)

    del exp_train_data['Description']
    del ess_train_data['Description']
    del exp_board_data['Description']

    return exp_train_data, ess_train_data, exp_board_data