#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate my-rdkit-env

PYARGS=""
PYARGS="$PYARGS --name Ep0_rewFacNone_stochastic_kernel"
PYARGS="$PYARGS --data_path src/dataset"
PYARGS="$PYARGS --artifact_path artifact"
PYARGS="$PYARGS --gpu 0"
PYARGS="$PYARGS --surrogate_model_url https://portal.nersc.gov/project/m3623/docking_score_models/NSP15_6W01_A_1_F_uncharged_upsamp/predict_logp/best_model.pth"
PYARGS="$PYARGS --surrogate_reward_timestep_delay 10"
PYARGS="$PYARGS --prob_redux_factor 0.99"
PYARGS="$PYARGS --stochastic_kernel"
PYARGS="$PYARGS --layer_num_g 4"

python src/main_gcpn.py $PYARGS
