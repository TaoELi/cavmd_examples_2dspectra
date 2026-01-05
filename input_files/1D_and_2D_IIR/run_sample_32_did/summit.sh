#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=qsummit
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --array=1-512
#SBATCH --output=./waste/nve_%A_%a.out
#SBATCH --error=./waste/nve_%A_%a.err

newE0=0
E0=$(echo "print(f'{($newE0)*(2**0.5):.4f}')" | python)
traj=$SLURM_ARRAY_TASK_ID
CURRENTFOLDER=$(pwd)
echo nve_"$traj"
cp $CURRENTFOLDER/data/nvt.sh ./waste/nvt_"$newE0"_"$traj".sh
#cp $CURRENTFOLDER/data/nve.sh ./waste/nvt_"$newE0"_"$traj".sh
sed -i "s/E0=0/E0=$E0/" ./waste/nvt_"$newE0"_"$traj".sh
sed -i "s/traj=0/traj=$traj/" ./waste/nvt_"$newE0"_"$traj".sh
sh ./waste/nvt_"$newE0"_"$traj".sh
