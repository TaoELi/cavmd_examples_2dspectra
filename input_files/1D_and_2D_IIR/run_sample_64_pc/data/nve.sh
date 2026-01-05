#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --job-name=nve
#SBATCH --partition=standard
#SBATCH --time=7-00:00:00
#SBATCH --output=./waste/%x-%j.out
#SBATCH --error=./waste/%x-%j.out

E0=0e-4
traj=0
echo "running equilibrium NVE for E0=$E0, traj=$traj"

HOMEFOLDER=$(pwd)
ORIGINFOLDER=$HOMEFOLDER/data/nve_sample
CHECKPOINTFOLDER=$HOMEFOLDER/data/init_water_xyz_nve

TOTALNVEFOLDER=$HOMEFOLDER/data/nve

if [ ! -d $TOTALNVEFOLDER ]; then
    mkdir $TOTALNVEFOLDER
fi

NVEFOLDER=$TOTALNVEFOLDER/nve_"$traj"

if [ ! -d $NVEFOLDER ]; then

    mkdir $NVEFOLDER
    cp $ORIGINFOLDER/* $NVEFOLDER
    sleep 1s
    cd $NVEFOLDER
    sed -i "s/water_280/h2o.E0_$E0.traj_$traj/" input.xml
    sed -i "s/water_280/h2o.E0_$E0.traj_$traj/" in.lmp
    sed -i "s/h2o-dipole/h2o-dipole_$E0.traj_$traj/" input.xml
    sed -i "s/RESTART/init_nve_$traj.checkpoint/" input.xml
    sed -i "s/'simu'/'simu_nve_$traj'/g" input.xml
    sed -i "s/0e-4/$E0/g" input.xml
    sed -i "s/<seed>31258<\/seed>/<seed>$RANDOM<\/seed>/" input.xml

fi

NEWNVECHECKPOINTFILE=$CHECKPOINTFOLDER/init_nve_$traj.checkpoint
cp $NEWNVECHECKPOINTFILE $NVEFOLDER
cd $NVEFOLDER

NVEPATTERN="<step>100000</step>"
CHECKFILE=$NVEFOLDER/simu_nve_"$traj".checkpoint

if [ -f $CHECKFILE ] && grep -q $NVEPATTERN "$CHECKFILE"; then
    echo nve_"$traj" finished
else
    echo nve_"$traj" not finished
    rm -rf /tmp/ipi_h2o.E0_$E0.traj_$traj
    rm -rf /tmp/ipi_h2o-dipole_$E0.traj_$traj
    if [ ! -f $CHECKFILE ] ; then
        i-pi input.xml & > output.output &
        echo "run new simulation"
    else
        i-pi $CHECKFILE & > output.output &
        echo "run existed simulation"
    fi
    sleep 10s
    lmp < in.lmp &> out_lammps &
    i-pi-driver -u -h h2o-dipole_$E0.traj_$traj -m water_dip_pol -o 1 &> out
    sleep 10s
fi

sed -i "s|<checkpoint stride='4' filename='checkpoint' overwrite='true'/>| |" $NVEFOLDER/input.xml
rm -rf /tmp/ipi_h2o.E0_$E0.traj_$traj
rm -rf /tmp/ipi_h2o-dipole_$E0.traj_$traj
#Run noneqm-traj script. Similar to i-pi but requires two inputs: i-pi input and spectra-related input.
./noneqm-traj.py input.xml -e 0.1 -t 1000 &> output &
sleep 10
lmp < in.lmp &> out_lammps &
i-pi-driver -u -h h2o-dipole_$E0.traj_$traj -m water_dip_pol -o 0 &> out
sleep 10

NONEQIRFOLDER=$HOMEFOLDER/result/noneq-ir-raman
if [ ! -d $NONEQIRFOLDER ]; then
    mkdir $NONEQIRFOLDER
fi

python ./noneqm-response.py input.xml -e 0.1 -t 1000 -b dip,pol -p 0,0
mv $NVEFOLDER/dip_pol_neq_2d.dat $NONEQIRFOLDER/dip_pol_neq_2d.dat_$traj
sleep 1s

DIPNVEFOLDER=$HOMEFOLDER/result/eqdip
if [ ! -d $DIPNVEFOLDER ]; then
    mkdir $DIPNVEFOLDER
fi

POLNVEFOLDER=$HOMEFOLDER/result/eqpol
if [ ! -d $POLNVEFOLDER ]; then
    mkdir $POLNVEFOLDER
fi

XYZNVEFOLDER=$HOMEFOLDER/result/eqxyz
if [ ! -d $XYZNVEFOLDER ]; then
    mkdir $XYZNVEFOLDER
fi

mv $NVEFOLDER/simu_nve_$traj.xc.xyz $XYZNVEFOLDER
mv $NVEFOLDER/simu_nve_$traj.dip_0 $DIPNVEFOLDER
mv $NVEFOLDER/simu_nve_$traj.pol_0 $POLNVEFOLDER