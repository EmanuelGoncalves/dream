__author__ = 'daniel'

from pandas import read_csv, DataFrame
from numpy import mean, std
from numpy.random import randn
from scipy.stats import spearmanr
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

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


def load_datasets():
    exp_train_data = load_gct_data(EXP_TRAINING_DATA)
    ess_train_data = load_gct_data(ESS_TRAINING_DATA)
    exp_board_data = load_gct_data(EXP_LEADERBOARD_DATA)

    del exp_train_data['Description']
    del ess_train_data['Description']
    del exp_board_data['Description']

    return exp_train_data, ess_train_data, exp_board_data


def random_essentiality():
    _, ess_train_data, _ = load_datasets()

    cell_lines = load_cell_lines(CELL_LINES_LEADERBOARD)
    null_model = lambda x: (mean(x) + randn(cell_lines.index.size) * std(x)).tolist()
    rand_data = [null_model(row) for _, row in ess_train_data.iterrows()]

    ess_rand_data = DataFrame(rand_data,
                              index=ess_train_data.index,
                              columns=cell_lines.index)

    ess_rand_data.insert(0, 'Description', ess_train_data.index)

    save_gct_data(ess_rand_data, 'rand_essentiality.gct')


def average_by_cell_line():
    _, ess_train_data, _ = load_datasets()

    lines_board = load_cell_lines(CELL_LINES_LEADERBOARD)
    lines_train = load_cell_lines(CELL_LINES_TRAINING)

    data = {}

    for line in lines_board.index:
        site = lines_board.at[line, 'Site_primary']
        matches = lines_train.index[lines_train['Site_primary'] == site]
        if matches.size > 0:
            data[line] = ess_train_data.loc[:, matches].mean(1).tolist()
        else:
            data[line] = ess_train_data.mean(1).tolist()

    ess_avg_data = DataFrame(data,
                             index=ess_train_data.index,
                             columns=lines_board.index)

    ess_avg_data.insert(0, 'Description', ess_train_data.index)
    save_gct_data(ess_avg_data, 'avg_per_line.gct')



def pre_process_data(X, Z, method='id', method_args={}):
    X2, Z2 = [], []

    if method == 'id':
        X2, Z2 = X, Z

    if method == 'z-score':
        z_score = X.std(0) / X.mean(0)
        X2 = X[:, z_score > method_args['z_min']]
        Z2 = Z[:, z_score > method_args['z_min']]


    scaler = StandardScaler().fit(X2)
    X2 = scaler.transform(X2)
    Z2 = scaler.transform(Z2)


    return X2, Z2


def train_and_predict(X, Y, Z, method, method_args):
    W = []
    Y2 = []

    if method == 'knn':
        estimator = KNeighborsRegressor(**method_args)

    if method == 'svm':
        estimator = SVR(**method_args)

    print '------------------>'
    for i, y in enumerate(Y):
        estimator.fit(X, y)
        y2 = estimator.predict(X)
        w = estimator.predict(Z)
        W.append(w)
        Y2.append(y2)
        if (i+1) % (len(Y) / 10) == 0:
            print '.',
    print ''


    return W, Y2


def training_score(Y, Y2):
    return mean([spearmanr(y, y2)[0] for y, y2 in zip(Y, Y2)])


def run_pipeline(preprocess, method, outputfile, pre_process_args={}, method_args={}):
    exp_train_data, ess_train_data, exp_board_data = load_datasets()

    #ess_train_data = ess_train_data.head(100)

    X = exp_train_data.values.T
    Y = ess_train_data.values
    Z = exp_board_data.values.T

    print 'pre-processing data using method', preprocess
    X2, Z2 = pre_process_data(X, Z, preprocess, pre_process_args)

    print 'training and predicting using method', method
    W, Y2 = train_and_predict(X2, Y, Z2, method, method_args)

    score = training_score(Y, Y2)
    print 'score for training dataset:', score

    ess_board_data = DataFrame(W,
                               columns=exp_board_data.columns,
                               index=ess_train_data.index)

    ess_board_data.insert(0, 'Description', ess_train_data.index)
    save_gct_data(ess_board_data, outputfile)

def main():
    #random_essentiality()
    #average_by_cell_line()
    #run_pipeline('z-score', 'knn', 'z_min_02_knn_5.gct', {'z_min': 0.2}, {'n_neighbors': 5})
    run_pipeline('z-score', 'svm', 'z_min_005_svm.gct', {'z_min': 0.05}, {})


if __name__ == '__main__':
    main()