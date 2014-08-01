__author__ = 'daniel'

from pandas import read_csv, DataFrame
from numpy import mean, std
from numpy.random import randn

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
            file2.write('#1.0\n')
            file2.write('{}\t{}\n'.format(dataframe.shape[0], dataframe.shape[1] - 1))
            for line in file1:
                file2.write(line)


def load_cell_lines(filename):
    data = read_csv(filename, sep='\t', header=0, index_col=0)
    return data


def random_essentiality():
    ess_train_data = load_gct_data(ESS_TRAINING_DATA)
    del ess_train_data['Description']

    cell_lines = load_cell_lines(CELL_LINES_LEADERBOARD)
    null_model = lambda x: (mean(x) + randn(cell_lines.index.size) * std(x)).tolist()
    rand_data = [null_model(row) for _, row in ess_train_data.iterrows()]

    ess_rand_data = DataFrame(rand_data,
                             index=ess_train_data.index,
                             columns=cell_lines.index)

    ess_rand_data.insert(0, 'Description', ess_train_data.index)

    save_gct_data(ess_rand_data, 'rand_essentiality.gct')


def average_by_cell_line():
    ess_train_data = load_gct_data(ESS_TRAINING_DATA)
    del ess_train_data['Description']

    lines_board = load_cell_lines(CELL_LINES_LEADERBOARD)
    lines_train = load_cell_lines(CELL_LINES_TRAINING)

    data = {}

    for line in lines_board.index:
        site = lines_board.at[line, 'Site_primary']
        matches = lines_train.index[lines_train['Site_primary'] == site]
        if matches.size > 0:
            data[line] = ess_train_data.loc[:, matches].apply(mean, 1).tolist()
        else:
            data[line] = ess_train_data.apply(mean, 1).tolist()

    ess_avg_data = DataFrame(data,
                             index=ess_train_data.index,
                             columns=lines_board.index)

    ess_avg_data.insert(0, 'Description', ess_train_data.index)

    save_gct_data(ess_avg_data, 'avg_per_line.gct')


def main():
    random_essentiality()
    average_by_cell_line()


if __name__ == '__main__':
    main()