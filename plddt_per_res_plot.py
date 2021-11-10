# Description: Plotting AlphaFold pLDDT score per residue
# Author: Vivian Monzon
# Adapted from sokrypton, Colabfold, 2021, https://github.com/sokrypton/ColabFold
# Mirdita M, Ovchinnikov S, Steinegger M. ColabFold - Making protein folding accessible to all. bioRxiv, 2021

import matplotlib.pyplot as plt
import pickle
import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os.path

parser = argparse.ArgumentParser()
parser.add_argument('--model_num', default=5, help='number of models')
parser.add_argument('--pkl_file_dir', help='path to alphafold results folder,'
                    'which contains the pickle ".pkl" files for each model.')
parser.add_argument('--query_seq', help='path to query sequence')
args = parser.parse_args()


def retieve_plddt_values(model_name):
    with open('{}/result_model_{}_ptm.pkl'.format(args.pkl_file_dir, model_name), 'rb') as f:
        data = pickle.load(f)
        plddt_list = data['plddt'].tolist()
    return plddt_list


model_plddts = {}
for x in range(1, int(args.model_num) + 1):
    if os.path.isfile('{}/result_model_{}_ptm.pkl'.format(
            args.pkl_file_dir, x)):
        plddt_per_model = retieve_plddt_values(x)
        model_plddts[x] = plddt_per_model
    else:
        continue


fh_query_seq = open(args.query_seq)
for name, seq in SimpleFastaParser(fh_query_seq):
    query_sequence = seq
    jobname = name.split(' ')[0]
    
homooligomer = 1


plt.title("Predicted lDDT per position ({})".format(jobname))
for m in range(1, int(args.model_num) + 1):
    if os.path.isfile('{}/result_model_{}_ptm.pkl'.format(
            args.pkl_file_dir, m)):
        plt.plot(model_plddts[m], label='{}_model'.format(m))
    else:
        continue
if homooligomer > 0:
    for n in range(homooligomer+1):
        x = n*(len(query_sequence)-1)
        plt.plot([x, x], [0, 100], color="black")
plt.legend()
plt.ylim(0, 100)
plt.ylabel("Predicted lDDT")
plt.xlabel("Positions")
plt.savefig('{}/{}'.format(args.pkl_file_dir, jobname+"_ptm_plDDT.png"),
            dpi=600)
