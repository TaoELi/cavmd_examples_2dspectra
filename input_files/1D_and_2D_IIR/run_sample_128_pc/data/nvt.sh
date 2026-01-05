#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --job-name=nvt
#SBATCH --partition=standard
#SBATCH --time=7-00:00:00
#SBATCH --output=./waste/%x-%j.out
#SBATCH --error=./waste/%x-%j.out

E0=0e-4
traj=0
echo "running NVT for E0=$E0, traj=$traj"

HOMEFOLDER=$(pwd)
ORIGINFOLDER=$HOMEFOLDER/data/nvt_sample
CHECKPOINTFOLDER=$HOMEFOLDER/data/init_water_xyz_nvt
NEWCHECKPOINTFOLDER=$HOMEFOLDER/data/init_water_xyz_nve

if [ ! -d $NEWCHECKPOINTFOLDER ]; then
    mkdir $NEWCHECKPOINTFOLDER
fi

TOTALNVTFOLDER=$HOMEFOLDER/data/nvt

if [ ! -d $TOTALNVTFOLDER ]; then
    mkdir $TOTALNVTFOLDER
fi

# run 1 noneq NVT simulation
NVTFOLDER=$TOTALNVTFOLDER/nvt_"$traj"

if [ ! -d $NVTFOLDER ]; then

    mkdir $NVTFOLDER
    cp $ORIGINFOLDER/* $NVTFOLDER
    sleep 1s
    cd $NVTFOLDER
    sed -i "s/0e-4/$E0/g" input.xml
    sed -i "s/water_280/h2o.nvt.E0_"$E0".traj_"$traj"/" input.xml
    sed -i "s/water_280/h2o.nvt.E0_"$E0".traj_"$traj"/" in.lmp
    sed -i "s/water_64.xyz/water_64.xyz_"$traj"/" input.xml
    sed -i "s/'simu'/'simu_nvt_"$traj"'/g" input.xml
    sed -i "s/<seed>31258<\/seed>/<seed>$RANDOM<\/seed>/" input.xml

fi

NEWNVECHECKPOINTFILE=$CHECKPOINTFOLDER/water_64.xyz_$traj
cp $NEWNVECHECKPOINTFILE $NVTFOLDER
cd $NVTFOLDER

NVTPATTERN="<step>40000</step>"
CHECKFILE=$NVTFOLDER/simu.init_nve
if [ -f $CHECKFILE ] && grep -q $NVTPATTERN "$CHECKFILE"; then
    echo nvt_"$traj" finished
else
    echo nvt_"$traj" not finished
    rm -rf /tmp/ipi_h2o.nvt.E0_"$E0".traj_"$traj"
    if [ ! -f $CHECKFILE ] ; then
        i-pi input.xml & > output.output &
        echo "run new simulation"
    else
        i-pi $CHECKFILE & > output.output &
        echo "run existed simulation"
    fi
    sleep 10s
    lmp < in.lmp
    sleep 10s
        
    cp simu.init_nve init_nve_"$traj".checkpoint
    mv $NVTFOLDER/init_nve_"$traj".checkpoint $NEWCHECKPOINTFOLDER

fi
