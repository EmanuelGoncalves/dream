__author__ = 'daniel'

from zipfile import ZipFile
from pandas import read_csv
from synapseclient import Synapse, Project, File

PROJECT_ID = 'syn2549343'
TEAM = 'UM-EBI'

#CODES = {'sc1': '2468319', 'sc2': '2468322', 'sc3': '2482339'}
CODES = {'sc1': '2571160', 'sc2': '2571162', 'sc3': '2571164'}

RESULTS_FOLDER = '../submissions/'


PHASE1_DATA = '../data/phase2/'
PHASE2_DATA = '../data/phase2/'

CNV_TRAINING_DATA_PH1 = PHASE1_DATA + 'CCLE_copynumber_training.gct'
EXP_TRAINING_DATA_PH1 = PHASE1_DATA + 'CCLE_expression_training.gct'
ESS_TRAINING_DATA_PH1 = PHASE1_DATA + 'Achilles_v2.9_training.gct'
CELL_LINES_TRAINING_PH1 = PHASE1_DATA + 'Cell_line_annotation_training.txt'
CNV_LEADERBOARD_DATA_PH1 = PHASE1_DATA + 'CCLE_copynumber_leaderboard.gct'
EXP_LEADERBOARD_DATA_PH1 = PHASE1_DATA + 'CCLE_expression_leaderboard.gct'
CELL_LINES_LEADERBOARD_PH1 = PHASE1_DATA + 'Cell_line_annotation_leaderboard.txt'
PRIORITY_GENE_LIST_PH1 = PHASE1_DATA + 'prioritized_gene_list.txt'

CNV_TRAINING_DATA_PH2 = PHASE2_DATA + 'CCLE_copynumber_training_phase2.gct'
EXP_TRAINING_DATA_PH2 = PHASE2_DATA + 'CCLE_expression_training_phase2.gct'
ESS_TRAINING_DATA_PH2 = PHASE2_DATA + 'Achilles_v2.11_training_phase2.gct'
CNV_LEADERBOARD_DATA_PH2 = PHASE2_DATA + 'CCLE_copynumber_leaderboard_phase2.gct'
EXP_LEADERBOARD_DATA_PH2 = PHASE2_DATA + 'CCLE_expression_leaderboard_phase2.gct'

CELL_LINES_TRAINING = PHASE2_DATA + 'Cell_line_annotation_training_phase2.txt'
CELL_LINES_LEADERBOARD = PHASE2_DATA + 'Cell_line_annotation_leaderboard_phase2.txt'
PRIORITY_GENE_LIST = PHASE2_DATA + 'prioritized_gene_list_phase2.txt'


def load_gct_data(filename):
    data = read_csv(filename, sep='\t', header=2, index_col=0)
    del data['Description']
    return data


def save_gct_data(dataframe, filename, folder=RESULTS_FOLDER):
    dataframe.to_csv(folder + 'temp.csv', sep='\t')
    with open(folder + 'temp.csv', 'r') as file1:
        with open(folder + filename, 'w') as file2:
            file2.write('#1.0\n')
            file2.write('{0}\t{1}\n'.format(dataframe.shape[0], dataframe.shape[1] - 1))
            for line in file1:
                file2.write(line)


def load_cell_lines(filename):
    data = read_csv(filename, sep='\t', header=0, index_col=0)
    return data

def load_gene_list(filename):
    data = read_csv(filename, header=None)
    return data.values[:,0]


def load_datasets(phase='phase2'):

    if phase == 'phase1':
        exp_train_data = load_gct_data(EXP_TRAINING_DATA_PH1)
        ess_train_data = load_gct_data(ESS_TRAINING_DATA_PH1)
        exp_board_data = load_gct_data(EXP_LEADERBOARD_DATA_PH1)

    if phase == 'phase2':
        exp_train_data = load_gct_data(EXP_TRAINING_DATA_PH2)
        ess_train_data = load_gct_data(ESS_TRAINING_DATA_PH2)
        exp_board_data = load_gct_data(EXP_LEADERBOARD_DATA_PH2)

    return exp_train_data, ess_train_data, exp_board_data

def submit_to_challenge(filename, challenge, label):

    client = Synapse()
    client.login()
    evaluation = client.getEvaluation(CODES[challenge])
#    client.joinEvaluation(evaluation)
    filename = filename + '.gct' if challenge == 'sc1' else filename + '.zip'
    myfile = File(RESULTS_FOLDER + filename, parent=PROJECT_ID)
    myfile = client.store(myfile)
    client.submit(evaluation, myfile, name=label, teamName=TEAM)

def zip_files(output, filelist):
    zipfile = ZipFile(RESULTS_FOLDER + output + '.zip', 'w')
    for filename in filelist:
        zipfile.write(RESULTS_FOLDER + filename, arcname=filename)
    zipfile.close()
