#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --job-name=analysis
#SBATCH --partition=standard
#SBATCH --time=7-00:00:00
#SBATCH --output=./%x-%j.out
#SBATCH --error=./%x-%j.out

HOMEFOLDER=$(pwd)/result
AVERAGEFILE=$HOMEFOLDER/average.py
FINALFOLDER=$HOMEFOLDER/final_results

if [ ! -f $FINALFOLDER ] ; then
    mkdir $FINALFOLDER
fi

sed -i "s|input=None|input='$HOMEFOLDER'|" $AVERAGEFILE
sed -i "s|output=None|output='$FINALFOLDER'|" $AVERAGEFILE

python -u $AVERAGEFILE
