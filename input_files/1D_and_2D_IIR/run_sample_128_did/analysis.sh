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

XYZNVEFOLDER=$HOMEFOLDER/eqxyz
PCDACFOLDERDS=$HOMEFOLDER/dac_simple_charge
if [ ! -f $PCDACFOLDERDS ] ; then
    mkdir $PCDACFOLDERDS
fi
ANALYSISFILE=$HOMEFOLDER/dac_simple_charge.py
python -u $ANALYSISFILE $XYZNVEFOLDER
sleep 3s
mv $XYZNVEFOLDER/*dac* $PCDACFOLDERDS

PHPCDACFOLDERDS=$HOMEFOLDER/phdac_simple_charge
if [ ! -f $PHPCDACFOLDERDS ] ; then
    mkdir $PHPCDACFOLDERDS
fi
ANALYSISFILE=$HOMEFOLDER/phac_simple_charge.py
python -u $ANALYSISFILE $XYZNVEFOLDER
sleep 3s
mv $XYZNVEFOLDER/*phdac* $PHPCDACFOLDERDS

NEWXYZNVEFOLDER=$HOMEFOLDER/eqdip
DACFOLDERDS=$HOMEFOLDER/dac_dipole_surface
if [ ! -f $DACFOLDERDS ] ; then
    mkdir $DACFOLDERDS
fi
NEWANALYSISFILE=$HOMEFOLDER/dac_dipole_surface.py
python -u $NEWANALYSISFILE $NEWXYZNVEFOLDER
sleep 3s
mv $NEWXYZNVEFOLDER/*dac* $DACFOLDERDS

NEWNEWXYZNVEFOLDER=$HOMEFOLDER/eqpol
PACFOLDERDS=$HOMEFOLDER/pac_iso_dipole_surface
if [ ! -f $PACFOLDERDS ] ; then
    mkdir $PACFOLDERDS
fi
NEWNEWANALYSISFILE=$HOMEFOLDER/pac_iso_dipole_surface.py
python -u $NEWNEWANALYSISFILE $NEWNEWXYZNVEFOLDER
sleep 3s
mv $NEWNEWXYZNVEFOLDER/*pac* $PACFOLDERDS

ANISOPACFOLDERDS=$HOMEFOLDER/pac_aniso_dipole_surface
if [ ! -f $ANISOPACFOLDERDS ] ; then
    mkdir $ANISOPACFOLDERDS
fi
ANISOANALYSISFILE=$HOMEFOLDER/pac_aniso_dipole_surface.py
python -u $ANISOANALYSISFILE $NEWNEWXYZNVEFOLDER
sleep 3s
mv $NEWNEWXYZNVEFOLDER/*pac* $ANISOPACFOLDERDS

AVERAGEFILE=$HOMEFOLDER/average.py
FINALFOLDER=$HOMEFOLDER/final_results

if [ ! -f $FINALFOLDER ] ; then
    mkdir $FINALFOLDER
fi

sed -i "s|input=None|input='$HOMEFOLDER'|" $AVERAGEFILE
sed -i "s|output=None|output='$FINALFOLDER'|" $AVERAGEFILE

python -u $AVERAGEFILE
