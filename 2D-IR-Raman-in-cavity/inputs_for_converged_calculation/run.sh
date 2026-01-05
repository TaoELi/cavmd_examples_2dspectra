#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --job-name=nve
#SBATCH --partition=standard
#SBATCH --time=7-00:00:00
#SBATCH --output=./%x-%j.out
#SBATCH --error=./%x-%j.err

E0=7e-4
traj=1
HOMEFOLDER=$(pwd)
rm -rf /tmp/ipi_h2o.E0_$E0.traj_$traj
rm -rf /tmp/ipi_h2o-dipole_$E0.traj_$traj
i-pi input.xml & > output.output &
sleep 10s
lmp < in.lmp &> out_lammps &
i-pi-driver -u -h h2o-dipole_$E0.traj_$traj -m water_dip_pol -o 1 &> out
sleep 10s

rm -rf /tmp/ipi_h2o.E0_$E0.traj_$traj
rm -rf /tmp/ipi_h2o-dipole_$E0.traj_$traj
#Run noneqm-traj script. Similar to i-pi but requires two inputs: i-pi input and spectra-related input.
./noneqm-traj.py input.xml -e 0.1 -t 1000 &> output &
sleep 10
lmp < in.lmp &> out_lammps &
i-pi-driver -u -h h2o-dipole_$E0.traj_$traj -m water_dip_pol -o 0 &> out
sleep 10

python ./noneqm-response.py input.xml -e 0.1 -t 1000 -b dip,pol -p 0,0
