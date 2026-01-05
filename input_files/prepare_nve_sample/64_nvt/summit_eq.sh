#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=qsummit
#SBATCH --partition=standard
#SBATCH --time=7-00:00:00
#SBATCH --array=1-10
#SBATCH --output=./nve_%A_%a.out
#SBATCH --error=./nve_%A_%a.err

traj=$SLURM_ARRAY_TASK_ID
CURRENTFOLDER=$(pwd)
echo nve_"$traj"
cp $CURRENTFOLDER/64_h2o/eqh2o.sh $CURRENTFOLDER/eqh2o_"$traj".sh
sed -i "s/nmol=100/nmol=$traj/" $CURRENTFOLDER/eqh2o_"$traj".sh
sh $CURRENTFOLDER/eqh2o_"$traj".sh
