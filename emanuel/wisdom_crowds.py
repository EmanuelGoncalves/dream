__author__ = 'emanuel'

from elastic_net import read_gct, save_gct_data

# Folders
submissions_folder = 'submissions/'

submissions_filenames = ['umebi_emanuel_13.gct',
                         'umebi_emanuel_14.gct',
                         'umebi_emanuel_15.gct',
                         'umebi_emanuel_16.gct',
                         'umebi_emanuel_17.gct',
                         'umebi_emanuel_18.gct',
                         'umebi_emanuel_19.gct',
                         'umebi_emanuel_20.gct']

submissions_dataframes = []

for filename in submissions_dataframes:
    read_gct(submissions_folder + filename)