#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate gcpndock

DATA=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

PYARGS=""
PYARGS="$PYARGS --name greedy_eval"
PYARGS="$PYARGS --data_path $DATA/src/dataset"
PYARGS="$PYARGS --warm_start_dataset_path /global/home/users/adchen/MD/2col/uncharged_unique/NSP15_6W01_A_1_F_combined_sorted_negonly_unique.csv"
PYARGS="$PYARGS --greedy"
PYARGS="$PYARGS --reward_type logp"
PYARGS="$PYARGS --artifact_path $DATA/artifact/gcpn"
PYARGS="$PYARGS --surrogate_guide_path $DATA/artifact/test_surrogate.pth"
PYARGS="$PYARGS --surrogate_eval_path $DATA/artifact/test_surrogate.pth"
PYARGS="$PYARGS --gcpn_path /global/home/users/adchen/exaLearnMol/artifact/gcpn/saves/crem_parallel_GPU_8_default_surr_2021.04.29_16:16:15/02112_gcpn.pth"
# PYARGS="$PYARGS --greedy"

python src/main_evaluate.py $PYARGS
