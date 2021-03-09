#!/bin/bash

if [ -z "$1" ]; then
  TIND=nan
else
  TIND=${1}
fi

python3 ~/simulation/diffractionSimulation/diffraction.py --run "sim_validation" --molName "NO2" --xyz_file "NO2.xyz" --basis_folder "/cds/group/ued/scratch/N2O/axis_distributions/NO2/A/temp-100K" --time_ind ${TIND}
