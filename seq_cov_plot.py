# Description: Plotting AlphaFold msa sequence coverrage per residue
# Author: Vivian Monzon
# Adapted from sokrypton, Colabfold, 2021, https://github.com/sokrypton/ColabFold
# Mirdita M, Ovchinnikov S, Steinegger M. ColabFold - Making protein folding accessible to all. bioRxiv, 2021

# Requirement: fasta2aln (https://github.com/yfukasawa/piaco/tree/master/PSICOV/fasta2aln)

import os
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--msa_dir', help='path to alignment files,'
                    'e.g. ANF00001/msas')
parser.add_argument('--fasta2aln_dir', default='~/Documents/Scripts/fasta2aln',
                    help='path to fasta2aln script '
                    '(https://github.com/yfukasawa/piaco/tree/master/PSICOV/fasta2aln)')
parser.add_argument('--query_seq', help='path to query sequence')
parser.add_argument('--out_dir', help='output folder')
args = parser.parse_args()

fh_query_seq = open(args.query_seq)
for name, seq in SimpleFastaParser(fh_query_seq):
    query_sequence = seq
    jobname = name.split(' ')[0]


def sto_to_dict(stockholm_fh):
    records = SeqIO.parse('{}/{}.sto'.format(args.msa_dir, stockholm_fh),
                          'stockholm')
    SeqIO.write(records, '{}/{}.fa'.format(args.msa_dir, stockholm_fh),
                'fasta')
    cmd = args.fasta2aln_dir + ' ' + args.msa_dir + '/' + stockholm_fh + \
        '.fa' + ' > ' + args.msa_dir + '/' + stockholm_fh + '_adapt.fa'
    os.system(cmd)
    with open('{}/{}_adapt.fa'.format(args.msa_dir, stockholm_fh)) as file_in:
        sequences = []
        for line in file_in:
            line = line.strip('\n')
            sequences.append(line)
    stockholm_format = open('{}/{}.sto'.format(args.msa_dir, stockholm_fh))
    seq_w_ids = {}
    sequence_ids = []
    for line in stockholm_format:
        if line.strip() and not line.startswith(('#', '//')):
            sequence_ids.append(line.split()[0])
    seq_w_ids = dict(zip(sequence_ids, sequences))
    return seq_w_ids


uniref_dict = sto_to_dict('uniref90_hits')
mgnify_dict = sto_to_dict('mgnify_hits')

cmd_bfd = args.fasta2aln_dir + ' ' + args.msa_dir + '/bfd_uniclust_hits.a3m' +\
    ' > ' + args.msa_dir + '/bfd_uniclust_hits_adapt.fa'
os.system(cmd_bfd)
bfd_headers = []
bfd_fh = open('{}/bfd_uniclust_hits.a3m'.format(args.msa_dir))
for name, seq in SimpleFastaParser(bfd_fh):
    bfd_headers.append(name.split(' ')[0])
with open('{}/bfd_uniclust_hits_adapt.fa'.format(args.msa_dir)) as file_in:
    bfd_seqs = []
    for line in file_in:
        bfd_seqs.append(line.strip('\n'))
bfd_seqs_w_ids = dict(zip(bfd_headers, bfd_seqs))

all_msa_seqs = {**uniref_dict, **mgnify_dict, **bfd_seqs_w_ids}
all_seqs_list = []
for key, value in all_msa_seqs.items():
    all_seqs_list.append(value)
msa_arr = np.array([list(seq) for seq in all_seqs_list])
seqid = (np.array(list(query_sequence)) == msa_arr).mean(-1)
seqid_sort = seqid.argsort()
non_gaps = (msa_arr != "-").astype(float)
non_gaps[non_gaps == 0] = np.nan

plt.title("Sequence coverage ({})".format(jobname))
plt.imshow(non_gaps[seqid_sort]*seqid[seqid_sort, None],
           interpolation='nearest', aspect='auto',
           cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
plt.plot((msa_arr != "-").sum(0), color='black')
plt.xlim(-0.5, msa_arr.shape[1]-0.5)
plt.ylim(-0.5, msa_arr.shape[0]-0.5)
plt.colorbar(label="Sequence identity to query",)
plt.xlabel("Positions")
plt.ylabel("Sequences")
plt.savefig("{}/{}_msa_coverage.png".format(args.out_dir, jobname), dpi=600)
