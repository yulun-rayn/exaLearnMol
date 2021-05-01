#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate gcpndock

DATA=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

PYARGS=""
PYARGS="$PYARGS --name gcpn_logp"
PYARGS="$PYARGS --gpu 2"
PYARGS="$PYARGS --nb_procs 8"
PYARGS="$PYARGS --data_path $DATA/src/dataset"
PYARGS="$PYARGS --reward_type logp"
PYARGS="$PYARGS --warm_start_dataset_path /clusterfs/csdata/data/MD/2col/uncharged_unique/NSP15_6W01_A_1_F_combined_sorted_negonly_unique.csv"
PYARGS="$PYARGS --artifact_path $DATA/artifact/gcpn"
PYARGS="$PYARGS --surrogate_model_path $DATA/artifact/test_surrogate.pth"

python src/main_gcpn.py $PYARGS
