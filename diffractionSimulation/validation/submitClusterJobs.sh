#!/bin/bash 

OUTDIR=./output/

if [ -z "$2" ]; then
  echo "ERROR SUBMITTING JOBS!!!   Must give the name of the molecule and largest L value to submit!"
  exit
fi

MOL=${1}
FILETORUN=run_${MOL}_simulation.sh
NJ=${2}

OUTPUTDIR=${OUTDIR}/logs/
for (( j=0; j<=$NJ; j++ ))
do
  for (( k=$(( -1*j )); k<=$j; k++ ))
  do
    lmk=$j",0,"$k
    echo "Submitting "${lmk}
    sbatch -p psanaq -o ${OUTPUTDIR}"job_"${lmk}".log" --wrap="python3 ../diffraction.py --run sim_validation --molecule CF2IBr --LMK $lmk --xyz_file CF2IBr.xyz --basis_folder /cds/group/ued/scratch/N2O/axis_distributions/NO2/A/temp-100K --time_ind 0 --cluster True"
    #sbatch -p psanaq --DefaultTime=INFINITE -t 10-10:00 -o ${OUTPUTDIR}"job_"${lmk}".log" --wrap="python3 ~/simulation/diffractionSimulation/diffraction.py --run sim_validation --molecule CF2IBr --LMK $lmk --xyz_file CF2IBr.xyz --basis_folder /cds/group/ued/scratch/N2O/axis_distributions/NO2/A/temp-100K --time_ind 0 --cluster True"
  done
done
