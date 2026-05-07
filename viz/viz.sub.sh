#!/bin/bash
#SBATCH --job-name=2d_fluid_viz
#SBATCH --output=JP-output.dat
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
## Propogate all environment vars to the job
#SBATCH --export=ALL

#conda init bash
. ~/.bashrc
conda info
conda activate DataAnalysis

cd $SLURM_SUBMIT_DIR/
# The next two lines are to provide more information about the environemnt
# that we are running in (and can be removed)
echo I am in $(pwd)
which python
python data_viewer_time.py

