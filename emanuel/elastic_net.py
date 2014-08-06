import numpy as np
from pandas import read_csv
from sklearn.metrics import r2_score
from sklearn.linear_model import ElasticNet

# Utilities function
def read_gct(filename):
    data = read_csv(filename, sep='\t', header=2, index_col=0)
    del data['Description']
    return data

# Data-sets files
train_exp_file = 'data/CCLE_expression_training.gct'
train_cnv_file = 'data/CCLE_copynumber_training.gct'
train_ess_file = 'data/Achilles_v2.9_training.gct'

leader_exp_file = 'data/CCLE_expression_leaderboard.gct'
leader_cnv_file = 'data/CCLE_copynumber_leaderboard.gct'

# Import data
train_exp = read_gct(train_exp_file)
train_cnv = read_gct(train_cnv_file)
train_ess = read_gct(train_ess_file)

leader_exp = read_gct(leader_exp_file)
leader_cnv = read_gct(leader_cnv_file)

# Perform elastic net
en = ElasticNet(alpha=0.1, l1_ratio=0.7)



