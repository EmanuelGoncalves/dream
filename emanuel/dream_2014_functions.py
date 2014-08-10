__author__ = 'emanuel'

import os
import numpy as np
from synapseclient import *
from pandas import read_csv, DataFrame

# Folders
submissions_folder = 'submissions/'


# Utilities function
def read_gct(filename):
    data = read_csv(filename, sep='\t', header=2, index_col=0)
    del data['Description']
    return data.T


def save_gct_data(dataframe, submission_filename_prefix, folder=submissions_folder, sep='\t'):
    n_submissions = len([x for x in os.listdir('submissions/') if x.startswith(submission_filename_prefix) and x.endswith('.gct')]) + 1
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


def write_features(gene_features, submission_filename_prefix, folder=submissions_folder, sep='\t'):
    n_submissions = len([x for x in os.listdir('submissions/') if x.startswith(submission_filename_prefix) and x.endswith('.txt')]) + 1
    filename = folder + submission_filename_prefix + str(n_submissions) + '.txt'

    with open(filename, 'w') as f:
        for k, v in gene_features.items():
            f.write(k)
            for feature in v:
                f.write(sep + feature)
            f.write('\n')

    return filename

def write_features_sc3(gene_features, submission_filename_prefix, folder=submissions_folder, sep='\t'):
    n_submissions = len([x for x in os.listdir('submissions/') if x.startswith(submission_filename_prefix) and x.endswith('.txt')]) + 1
    filename = folder + submission_filename_prefix + str(n_submissions) + '.txt'

    with open(filename, 'w') as f:
        f.write('\t'.join(gene_features))
        f.write('\n')

    return filename

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