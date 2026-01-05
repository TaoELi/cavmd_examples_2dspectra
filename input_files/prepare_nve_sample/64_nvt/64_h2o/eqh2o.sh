#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
##SBATCH --mem-per-cpu=1024M
#SBATCH --job-name=test
#SBATCH --partition=taoeli
##SBATCH --nodelist=r05n04
#SBATCH --time=7-00:00:00
#SBATCH --output=./waste/%x-%j.out
#SBATCH --error=./waste/%x-%j.out

nmol=100
echo running folder is h2o_"$nmol"

TOTALFOLDER=$(pwd)
if [ ! -d $TOTALFOLDER ]; then
    mkdir $TOTALFOLDER
fi

HOMEFOLDER=$TOTALFOLDER/h2o_"$nmol"
if [ ! -d $HOMEFOLDER ]; then
    mkdir $HOMEFOLDER
fi

ORIGINFOLDER=$TOTALFOLDER/64_h2o

# run NVT equilibrate
CHECKPOINTFOLDER=$HOMEFOLDER/checkpoint
if [ ! -d $CHECKPOINTFOLDER ]; then
    mkdir $CHECKPOINTFOLDER
fi

NVTPATTERN="<step>400000</step>"
NVTFOLDER=$HOMEFOLDER/nvt

if [ -f "$CHECKPOINTFOLDER/init_0.checkpoint" ]; then
    
    if grep -q $NVTPATTERN "$CHECKPOINTFOLDER/init_0.checkpoint"; then
        echo 'NVT equilibrate for 100 ps has finished'
    fi
    
else

    if [ ! -d "$NVTFOLDER" ]; then

        mkdir $NVTFOLDER
        cp $ORIGINFOLDER/* $NVTFOLDER
        cd $NVTFOLDER
        sed -i "s/<seed>31258<\/seed>/<seed>$RANDOM<\/seed>/" input_eq.xml
        sed -i "s/water_280/water_280_"$nmol"_nvt_64/" input_eq.xml
        sed -i "s/water_280/water_280_"$nmol"_nvt_64/" in.lmp

    fi
    
    cd $NVTFOLDER
    echo 'NVT checkpoint file does not exist, run NVT equilibrate 100 ps'
    rm -rf /tmp/ipi_water_280_"$nmol"_nvt_64
    i-pi input_eq.xml & > nvteq.output &
    sleep 5s
    lmp < in.lmp
    sleep 5s
            
    cd $HOMEFOLDER
    mv $NVTFOLDER/simu.init_nve $CHECKPOINTFOLDER/init_0.checkpoint

fi

# run 1 NVT simulation
NVEPATTERN="<step>256000</step>"
traj=1
    
    NVECHECKPOINTFILE=$CHECKPOINTFOLDER/init_$traj.checkpoint
    NVEFOLDER=$HOMEFOLDER/nvt_$traj
    
    if [ -f "$NVECHECKPOINTFILE" ]; then
    
        if grep -q $NVEPATTERN "$NVECHECKPOINTFILE"; then
    	      echo "Found checkpoint for $traj-th NVE simulation finished, skip $traj-th sequential job"
        fi

    else

        if [ ! -d "$NVEFOLDER" ]; then

            mkdir $NVEFOLDER
            cp $ORIGINFOLDER/* $NVEFOLDER
            cd $NVEFOLDER
            sed -i "s/<seed>31258<\/seed>/<seed>$RANDOM<\/seed>/" input_traj.xml
            sed -i "s/water_280/water_280_"$nmol"_nve_"$traj"_64/" input_traj.xml
            sed -i "s/water_280/water_280_"$nmol"_nve_"$traj"_64/" in.lmp
            sed -i "s/init_nve.checkpoint/init_$(($traj-1)).checkpoint/" input_traj.xml
                        
        fi
        
        NEWNVECHECKPOINTFILE=$CHECKPOINTFOLDER/init_$(($traj-1)).checkpoint
        cp $NEWNVECHECKPOINTFILE $NVEFOLDER
        cd $NVEFOLDER
        echo "$traj-th NVE checkpoint file does not exist, run NVE simulation"
        rm -rf /tmp/ipi_water_280_"$nmol"_nve_"$traj"_64
        i-pi input_traj.xml & > output.output &
        sleep 5s
        lmp < in.lmp
        sleep 5s
                        
    fi
    
