__author__ = 'emanuel'

from pandas import read_csv, DataFrame
import numpy as np
import synapseclient

# Folders
submissions_folder = 'submissions/'

def read_gct(filename):
    data = read_csv(filename, sep='\t', header=2, index_col=0)
    del data['Description']
    return data.T

def save_gct_data(dataframe, filename, folder=submissions_folder, sep='\t'):
    rows = dataframe.shape[0]
    cols = dataframe.shape[1]

    with open(folder + filename, 'w') as f:
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

submissions_filenames = ['umebi_emanuel_13.gct',
                         'umebi_emanuel_14.gct',
                         'umebi_emanuel_15.gct',
                         'umebi_emanuel_16.gct',
                         'umebi_emanuel_17.gct',
                         'umebi_emanuel_18.gct',
                         'umebi_emanuel_19.gct',
                         'umebi_emanuel_20.gct',
                         'umebi_emanuel_22.gct',
                         'umebi_emanuel_23.gct',
                         'umebi_emanuel_24.gct',
                         'umebi_emanuel_25.gct',
                         'umebi_emanuel_26.gct',
                         'umebi_emanuel_27.gct',
                         'umebi_emanuel_28.gct',
                         'umebi_emanuel_29.gct']

submissions_dataframes = []

for name in submissions_filenames:
    print name
    submissions_dataframes.append(read_gct(submissions_folder + name))

nrows, ncolumns = submissions_dataframes[0].shape
ndataframes = len(submissions_filenames)

averaged = DataFrame(None, index=submissions_dataframes[0].axes[1], columns=submissions_dataframes[0].axes[0])

for i in range(nrows):
    for j in range(ncolumns):
        averaged.ix[j, i] = np.mean([x.ix[i, j] for x in submissions_dataframes])

save_gct_data(averaged, 'umebi_emanuel_30.gct')

# Set-up synap client
syn = synapseclient.Synapse()
syn.login()
