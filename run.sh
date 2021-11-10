#!/bin/bash

main(){
    plddt_plots
    # seq_cov_plots
    plddt_seq_len
    pTM_seq_len
}

# antifam structure models provided under: https://drive.google.com/drive/folders/1u9OocRIAabGQn56GljoG1JTDAxjkY1ro
plddt_plots(){
    for i in $(seq -f "%03g" 1 251); do
    	python3.7 plddt_per_res_plot.py --model_num 5 --pkl_file_dir antifam_structure_models/ANF00${i}/ --query_seq query_seqs/ANF00${i}.fa
    done
}

# MSA folders not provided due to space requirements
seq_cov_plots(){
    for i in $(seq -f "%03g" 1 251); do
	python3.7 seq_cov_plot.py --msa_dir antifam/ANF00${i}/msas --query_seq query_seqs/ANF00${i}.fa --out_dir antifam_structure_models/ANF00${i}/
    done
}

plddt_seq_len(){
    python3.7 plddt_scores.py 
}

pTM_seq_len(){
    python3.7 pTM_scores.py
}

main
