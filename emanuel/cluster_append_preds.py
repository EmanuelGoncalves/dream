__author__ = 'emanuel'

import numpy as np
from pandas import DataFrame, Series
from dream_2014_functions import read_data_sets, save_gct_data, submit_solution, ev_code_sc1

# Folders
submission_filename_prefix = 'sc1_emanuel_phase2_'

# Import data-sets
train_exp, train_cnv, train_ess, leader_exp, leader_cnv, prioritized_genes = read_data_sets()

genes = train_ess.axes[1]
samples = leader_exp.axes[0]
predictions = DataFrame(None, index=genes, columns=samples)

for gene in genes:
    with open('_submissions_temp/' + gene + '.txt', 'r') as f:
        values = f.readline().split('\t')
        y_pred = Series(values).ix[range(2, len(values))]
        predictions.ix[values[0]] = np.array(y_pred, np.float64)

filename = save_gct_data(predictions, submission_filename_prefix)
print '[DONE]: Saved to file ' + filename

submit_solution(filename, filename.split('/')[1], ev_code_sc1)
print '[SUBMITED]'