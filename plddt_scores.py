# Description: Collecting pLDDT scores of Antifam structure predictions
# Author: Vivian Monzon
# Note: antifam_structure_models provided under: https://drive.google.com/drive/folders/1u9OocRIAabGQn56GljoG1JTDAxjkY1ro

import pickle
import statistics
from pathlib import Path
import pandas as pd


def retieve_plddt_mean(seq_name):
    with open('antifam_structure_models/{}/result_model_1_ptm.pkl'.format(
            seq_name), 'rb') as f:
        data = pickle.load(f)
        plddt_list = data['plddt'].tolist()
        plddt_mean = statistics.mean(plddt_list)
        seq_len = len(plddt_list)
    return plddt_mean, seq_len


name_seq = {}
name_plddt = {}
a = [num for num in range(1, 252)]
for x in a:
    name = 'ANF00' + str(x).rjust(3, '0')
    pkl_file = Path(
        'antifam_structure_models/{}/result_model_1_ptm.pkl'.format(name))
    if pkl_file.is_file():
        plddt_mean_value, seq_length = retieve_plddt_mean(name)
        name_seq[name] = seq_length
        name_plddt[name] = plddt_mean_value
    else:
        print('{} doesnot exists'.format(name))


df_plddt_mean = pd.DataFrame(name_plddt.items(),
                             columns=['Name', 'plDDT mean'])
df_seq_len = pd.DataFrame(name_seq.items(),
                          columns=['Name', 'Sequence length'])
df_merged = df_plddt_mean.merge(df_seq_len, on='Name', how='inner')
df_merged.to_csv('data/plddt_seq_length_info_ptm.csv', index=False)
