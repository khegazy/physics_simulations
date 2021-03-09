#!/bin/bash

if [ -z "$1" ]; then
  TIND=nan
else
  TIND=${1}
fi

#python3 ~/simulation/diffractionSimulation/diffraction.py --run "sim_validation" --molecule "CF2IBr" --xyz_file "CF2IBr.xyz" --basis_folder "/cds/group/ued/scratch/N2O/axis_distributions/NO2/A/temp-100K" --time_ind ${TIND} --cluster ${CLUSTER}
python3 ../diffraction.py --run "sim_validation" --molecule "CF2IBr" --LMK "20,0,1" --xyz_file "CF2IBr.xyz" --basis_folder "/cds/group/ued/scratch/N2O/axis_distributions/NO2/A/temp-100K" --time_ind ${TIND}
