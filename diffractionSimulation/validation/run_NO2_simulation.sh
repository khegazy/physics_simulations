#!/bin/bash

if [ -z "$1" ]; then
  TIME=nan
else
  TIME=${1}
fi

python3 ../diffraction.py --run "sim_validation" --molecule "NO2" --xyz_file "NO2.xyz" --basis_folder "/cds/group/ued/scratch/N2O/axis_distributions/NO2/A/temp-100K" --eval_time ${TIME}
