__author__ = 'daniel'

from zipfile import ZipFile
from pandas import read_csv
from synapseclient import Synapse, File

PROJECT_ID = 'syn2549343'
TEAM = 'UM-EBI'

#CODES = {'sc1': '2468319', 'sc2': '2468322', 'sc3': '2482339'}
#CODES = {'sc1': '2571160', 'sc2': '2571162', 'sc3': '2571164'}
CODES = {'sc1': '2571166', 'sc2': '2571168', 'sc3': '2571170'}

RESULTS_FOLDER = '../submissions/'

PHASE1_DATA = '../data/phase1/'
PHASE2_DATA = '../data/phase2/'
PHASE3_DATA = '../data/phase3/'

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
CELL_LINES_TRAINING_PH2 = PHASE2_DATA + 'Cell_line_annotation_training_phase2.txt'
CELL_LINES_LEADERBOARD_PH2 = PHASE2_DATA + 'Cell_line_annotation_leaderboard_phase2.txt'
PRIORITY_GENE_LIST_PH2 = PHASE2_DATA + 'prioritized_gene_list_phase2.txt'

CNV_TRAINING_DATA_PH3 = PHASE3_DATA + 'CCLE_copynumber_training_phase3.gct'
EXP_TRAINING_DATA_PH3 = PHASE3_DATA + 'CCLE_expression_training_phase3.gct'
ESS_TRAINING_DATA_PH3 = PHASE3_DATA + 'Achilles_v2.11_training_phase3.gct'
MUT_TRAINING_DATA_PH3 = PHASE3_DATA + 'CCLE_hybridmutation_training_phase3.gct'
CNV_LEADERBOARD_DATA_PH3 = PHASE3_DATA + 'CCLE_copynumber_finaltest_phase3.gct'
EXP_LEADERBOARD_DATA_PH3 = PHASE3_DATA + 'CCLE_expression_finaltest_phase3.gct'
MUT_LEADERBOARD_DATA_PH3 = PHASE3_DATA + 'CCLE_hybridmutation_finaltest_phase3.gct'
CELL_LINES_TRAINING_PH3 = PHASE3_DATA + 'Cell_line_annotation_training_phase3.txt'
CELL_LINES_LEADERBOARD_PH3 = PHASE3_DATA + 'Cell_line_annotation_finaltest_phase3.txt'
PRIORITY_GENE_LIST_PH3 = PHASE3_DATA + 'prioritized_gene_list_phase3.txt'


def load_gct_data(filename, index_by_description=False):
    if index_by_description:
        data = read_csv(filename, sep='\t', header=2, index_col=1)
        data = data.loc[[str(x) != 'nan' for x in data.index], :]
        del data['Name']
    else:
        data = read_csv(filename, sep='\t', header=2, index_col=0)
        del data['Description']
    return data


def save_gct_data(dataframe, filename, folder=RESULTS_FOLDER):
    dataframe.to_csv(folder + 'temp_' + filename, sep='\t')
    with open(folder + 'temp_' + filename, 'r') as file1:
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


def load_datasets(phase='phase2', get_cnv=False, get_mut=False):

    datasets = {}

    if phase == 'phase1':
        exp_train_data = EXP_TRAINING_DATA_PH1
        ess_train_data = ESS_TRAINING_DATA_PH1
        exp_board_data = EXP_LEADERBOARD_DATA_PH1
        cnv_train_data = CNV_TRAINING_DATA_PH1
        cnv_board_data = CNV_LEADERBOARD_DATA_PH1
        gene_list = PRIORITY_GENE_LIST_PH1

    if phase == 'phase2':
        exp_train_data = EXP_TRAINING_DATA_PH2
        ess_train_data = ESS_TRAINING_DATA_PH2
        exp_board_data = EXP_LEADERBOARD_DATA_PH2
        cnv_train_data = CNV_TRAINING_DATA_PH2
        cnv_board_data = CNV_LEADERBOARD_DATA_PH2
        gene_list = PRIORITY_GENE_LIST_PH2

    if phase == 'phase3':
        exp_train_data = EXP_TRAINING_DATA_PH3
        ess_train_data = ESS_TRAINING_DATA_PH3
        exp_board_data = EXP_LEADERBOARD_DATA_PH3
        cnv_train_data = CNV_TRAINING_DATA_PH3
        cnv_board_data = CNV_LEADERBOARD_DATA_PH3
        mut_train_data = MUT_TRAINING_DATA_PH3
        mut_board_data = MUT_LEADERBOARD_DATA_PH3
        gene_list = PRIORITY_GENE_LIST_PH3

    datasets['ess_train_data'] = load_gct_data(ess_train_data)

    datasets['exp_train_data'] = load_gct_data(exp_train_data)
    datasets['exp_board_data'] = load_gct_data(exp_board_data)

    # datasets['exp_train_data'] = load_gct_data(exp_train_data, True)
    # datasets['exp_board_data'] = load_gct_data(exp_board_data, True)
    # common_genes = set(datasets['exp_train_data'].index) & set(datasets['ess_train_data'].index)
    # datasets['exp_train_data'] = datasets['exp_train_data'].loc[common_genes, :]
    # datasets['exp_board_data'] = datasets['exp_board_data'].loc[common_genes, :]

    datasets['gene_list'] = load_gene_list(gene_list)

    if get_cnv:
        datasets['cnv_train_data'] = load_gct_data(cnv_train_data)
        datasets['cnv_board_data'] = load_gct_data(cnv_board_data)

    if phase == 'phase3' and get_mut:
        datasets['mut_train_data'] = load_gct_data(mut_train_data).fillna(0.5)
        datasets['mut_board_data'] = load_gct_data(mut_board_data).fillna(0.5)

    return datasets

def submit_to_challenge(filename, challenge, label, retry=True):

    try:
        client = Synapse()
        client.login()
        evaluation = client.getEvaluation(CODES[challenge])
        filename = filename + '.gct' if challenge == 'sc1' else filename + '.zip'
        myfile = File(RESULTS_FOLDER + filename, parent=PROJECT_ID)
        myfile = client.store(myfile)
        client.submit(evaluation, myfile, name=label, teamName=TEAM)
    except:
        if retry:
            submit_to_challenge(filename, challenge, label, retry=False)
        else:
            print 'failed to submit', label, 'to', challenge

def zip_files(output, filelist):
    zipfile = ZipFile(RESULTS_FOLDER + output + '.zip', 'w')
    for filename in filelist:
        zipfile.write(RESULTS_FOLDER + filename, arcname=filename)
    zipfile.close()
