# Description: Collecting pTM scores of Antifam structure predictions
# Author: Vivian Monzon
# Note: antifam_structure_models provided under: https://drive.google.com/drive/folders/1u9OocRIAabGQn56GljoG1JTDAxjkY1ro

import pickle
import pandas as pd
from pathlib import Path


def retieve_pTM(seq_name):
    with open('antifam_structure_models/{}/result_model_1_ptm.pkl'.format(
            seq_name), 'rb') as f:
        data = pickle.load(f)
        ptm = data['ptm']
        plddt_list = data['plddt'].tolist()
        seq_len = len(plddt_list)
    return ptm, seq_len


name_seq = {}
name_pTM = {}
a = [num for num in range(1, 252)]
for x in a:
    name = 'ANF00' + str(x).rjust(3, '0')
    pkl_file = Path(
        'antifam_structure_models/{}/result_model_1_ptm.pkl'.format(name))
    if pkl_file.is_file():
        pTM_value, seq_length = retieve_pTM(name)
        name_seq[name] = seq_length
        name_pTM[name] = pTM_value
    else:
        print('{} doesnot exists'.format(name))


df_pTM = pd.DataFrame(name_pTM.items(), columns=['Name', 'pTM'])
df_seq_len = pd.DataFrame(name_seq.items(),
                          columns=['Name', 'Sequence length'])
df_merged = df_pTM.merge(df_seq_len, on='Name', how='inner')

df_merged.to_csv('data/pTM_seq_length_info.csv', index=False)
