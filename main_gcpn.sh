#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate gcpndock

DATA=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

PYARGS=""
PYARGS="$PYARGS --name NSP15_6W01_A1F_crem_dock_2"
PYARGS="$PYARGS --data_path $DATA"
PYARGS="$PYARGS --artifact_path /clusterfs/csdata/data/gcpn"
PYARGS="$PYARGS --use_surrogate"
PYARGS="$PYARGS --use_crem"
PYARGS="$PYARGS --gpu 2"
PYARGS="$PYARGS --surrogate_model_1_path /clusterfs/csdata/data/surrogate_models/NSP15_combined1_upsample/predict_logp/best_model.pth"
PYARGS="$PYARGS --surrogate_model_2_path /clusterfs/csdata/data/surrogate_models/NSP15_combined2_upsample/predict_logp/best_model.pth"
PYARGS="$PYARGS --is_conditional 1"
PYARGS="$PYARGS --conditional /clusterfs/csdata/data/MD/2col/uncharged_unique/NSP15_6W01_A_1_F_combined_sorted_negonly_unique.csv"
PYARGS="$PYARGS --surrogate_reward_timestep_delay 10"

python3 src/main_gcpn.py $PYARGS
