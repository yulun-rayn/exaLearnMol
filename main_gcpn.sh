#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate my-rdkit-env_111_171_163

PYARGS=""
PYARGS="$PYARGS --name Ep0_rewFacNone_stochastic_kernel"
PYARGS="$PYARGS --data_path src/dataset"
PYARGS="$PYARGS --artifact_path aritfact"
PYARGS="$PYARGS --gpu 0"
PYARGS="$PYARGS --surrogate_model_path artifact/surrogate_model.pth"
PYARGS="$PYARGS --surrogate_reward_timestep_delay 10"
PYARGS="$PYARGS --prob_redux_factor 0.99"
PYARGS="$PYARGS --stochastic_kernel"
PYARGS="$PYARGS --layer_num_g 4"

python src/main_gcpn.py $PYARGS
