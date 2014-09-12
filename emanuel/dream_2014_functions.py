__author__ = 'emanuel'

import os
import numpy as np
from synapseclient import *
from pandas import read_csv

# Folders
submissions_folder = '_submissions_phase3/'

# Data-sets files
train_exp_file = '_data_phase3/CCLE_expression_training_phase3.gct'
train_cnv_file = '_data_phase3/CCLE_copynumber_training_phase3.gct'
train_ess_file = '_data_phase3/Achilles_v2.11_training_phase3.gct'

leader_exp_file = '_data_phase3/CCLE_expression_finaltest_phase3.gct'
leader_cnv_file = '_data_phase3/CCLE_copynumber_finaltest_phase3.gct'

prioritized_gene_list = '_data_phase3/prioritized_gene_list_phase3.txt'

# Annotation files
exp_annot_file = '_data_phase3/CCLE_expression_gene_annotations_phase3.txt'

# Evaluation codes
ev_code_sc1 = 2571160
ev_code_sc2 = 2571162
ev_code_sc3 = 2571164

# Cell lines phase 2
leader_board_cell_lines = [
    'CAL12T', 'CAL27', 'CAL54', 'COLO684', 'DB', 'DMS273', 'EBC1', 'ES2', 'G361', 'HCC15', 'HCC38', 'HSC3', 'KMS18',
    'LN18', 'MDAMB468', 'MEWO', 'NCIH1155', 'NCIH1355', 'NCIH1944', 'NCIH2291', 'NCIH460', 'NCIH522', 'NCIH526',
    'NCIH647', 'NCIH727', 'NUGC3', 'PC14', 'PC3', 'PECAPJ49', 'RD', 'SUDHL10', 'TUHR10TKB', 'UACC257']

training_cell_lines = ['143B', '769P', 'A3KAW', 'ABC1', 'CHAGOK1', 'CJM', 'COLO679', 'CORL105', 'DMS114', 'DND41',
    'G402', 'HCC1500', 'HCC1806', 'HCC366', 'HCC95', 'HEC251', 'HEC50B', 'HEKTE', 'HLC1', 'HUH1', 'INA6', 'JHH5',
    'JHOS2', 'KARPAS422', 'KMM1', 'KYSE70', 'LC1SQSF', 'LCLC103H', 'LCLC97TM1', 'LU99', 'M059K', 'MDAMB157', 'MFE280',
    'MKN45', 'NCIH1373', 'NCIH1395', 'NCIH146', 'NCIH1568', 'NCIH1581', 'NCIH1648', 'NCIH1703', 'NCIH1781', 'NCIH1836',
    'NCIH2009', 'NCIH2030', 'NCIH2073', 'NCIH209', 'NCIH2126', 'NCIH358', 'NCIH520', 'NCIH841', 'OCIMY7', 'OPM1',
    'OVTOKO', 'PLCPRF5', 'RH30', 'RL952', 'SCC4', 'SKMES1', 'SNU349', 'SNU878', 'SW1271', 'SW837', 'T47D', 'U2OS',
    'U937']


def read_gene_neighbours():
    genes_neighbours = {}
    with open('emanuel/list1_second_neighbours.txt', 'r') as f:
        lines = f.read().splitlines()

        for line in lines:
            gene = line.split('\t')[0]
            neighbours = line.split('\t')[1].split(';')

            genes_neighbours[gene] = neighbours

    return genes_neighbours


def read_annot(filename):
    data = read_csv(filename, sep='\t', index_col=0)
    return data


def read_annotations():
    return read_annot(exp_annot_file)


def read_gct(filename):
    data = read_csv(filename, sep='\t', header=2, index_col=0)
    del data['Description']
    return data.T


def read_data_sets():
    train_exp = read_gct(train_exp_file)
    train_cnv = read_gct(train_cnv_file)
    train_ess = read_gct(train_ess_file)

    leader_exp = read_gct(leader_exp_file)
    leader_cnv = read_gct(leader_cnv_file)

    prioritized_genes = np.genfromtxt(prioritized_gene_list, dtype='str')

    return train_exp, train_cnv, train_ess, leader_exp, leader_cnv, prioritized_genes


def save_gct_data(dataframe, submission_filename_prefix, folder=submissions_folder, sep='\t'):
    n_submissions = len([x for x in os.listdir(submissions_folder) if x.startswith(submission_filename_prefix) and x.endswith('.gct')]) + 1
    filename = folder + submission_filename_prefix + str(n_submissions) + '.gct'

    rows = dataframe.shape[0]
    cols = dataframe.shape[1]

    with open(filename, 'w') as f:
        f.write('#1.2\n')
        f.write(str(rows) + '\t' + str(cols) + '\n')

        f.write('Name' + sep + 'Description')
        for header in dataframe.axes[1]:
            f.write(sep + header)
        f.write('\n')

        for i in range(rows):
            row_name = dataframe.axes[0][i]
            f.write(row_name + sep + row_name)

            for j in range(cols):
                f.write(sep + '%.5f' % dataframe.ix[i, j])

            f.write('\n')

    return filename

def write_features(gene_features, submission_filename_prefix, folder=submissions_folder, sep='\t'):
    n_submissions = len([x for x in os.listdir(submissions_folder) if x.startswith(submission_filename_prefix) and x.endswith('.txt')]) + 1
    filename = folder + submission_filename_prefix + str(n_submissions) + '.txt'

    with open(filename, 'w') as f:
        for k, v in gene_features.items():
            f.write(k)
            for feature in v:
                f.write(sep + feature)
            f.write('\n')

    return filename


def write_features_sc3(gene_features, submission_filename_prefix, folder=submissions_folder, sep='\t'):
    n_submissions = len([x for x in os.listdir(submissions_folder) if x.startswith(submission_filename_prefix) and x.endswith('.txt')]) + 1
    filename = folder + submission_filename_prefix + str(n_submissions) + '.txt'

    with open(filename, 'w') as f:
        f.write(sep.join(gene_features))
        f.write('\n')

    return filename

def submit_solution(filepath, name, evaluation_id, team='UM-EBI', project_id='syn2563019'):
    syn = Synapse()
    syn.login()

    project = syn.get(project_id)
    entity = syn.store(File(filepath, parent=project))

    evaluation = syn.getEvaluation(evaluation_id)

    submission = syn.submit(evaluation, entity, name=name, teamName=team)

    print submission

    return submission


def resampling(y_train, X_train, nresampling=20):
    # Observations resampling
    q1 = y_train.quantile(0.25)
    q3 = y_train.quantile(0.75)
    iqr = q3 - q1

    is_outlier = np.logical_not(np.logical_and(y_train > (q1 - 1.5 * iqr), y_train < (q3 + 1.5 * iqr)))

    y_train_outliers = y_train[is_outlier]
    X_train_outliers = X_train[is_outlier.values]

    for i in range(nresampling):
        y_train = y_train.append(y_train_outliers)
        X_train = np.append(X_train, X_train_outliers, axis=0)


def read_features(filename):
    features_table = {}
    data = read_csv(filename, sep='\t', index_col=0, header=None)

    for gene in data.axes[0]:
        for feature in data.loc[gene]:
            if not features_table.has_key(feature):
                features_table[feature] = 1
            else:
                features_table[feature] += 1

    return features_table