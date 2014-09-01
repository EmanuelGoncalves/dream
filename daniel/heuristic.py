__author__ = 'daniel'

from pandas import DataFrame
from datasets import load_datasets, save_gct_data, load_cell_lines
from datasets import CELL_LINES_TRAINING_PH1, CELL_LINES_LEADERBOARD_PH1

def average_by_cell_line():
    datasets = load_datasets()
    ess_train_data = datasets['ess_train_data']

    lines_board = load_cell_lines(CELL_LINES_LEADERBOARD_PH1)
    lines_train = load_cell_lines(CELL_LINES_TRAINING_PH1)

    data = {}

    for line in lines_board.index:
        site = lines_board.at[line, 'Site_primary']
        matches = lines_train.index[lines_train['Site_primary'] == site]
        if matches.size > 0:
            data[line] = ess_train_data.loc[:, matches].mean(1).tolist()
        else:
            data[line] = ess_train_data.mean(1).tolist()

    ess_avg_data = DataFrame(data, index=ess_train_data.index, columns=lines_board.index)
    ess_avg_data.insert(0, 'Description', ess_train_data.index)
    save_gct_data(ess_avg_data, 'avg_per_line.gct')


