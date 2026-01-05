#!/bin/bash

E0=7e-4
traj=1
rm -rf /tmp/ipi_h2o.E0_$E0.traj_$traj
rm -rf /tmp/ipi_h2o-dipole_$E0.traj_$traj
i-pi input.xml & > output.output &
sleep 5
lmp < in.lmp &> out_lammps &
i-pi-driver -u -h h2o-dipole_$E0.traj_$traj -m water_dip_pol -o 1 &> out
sleep 5

rm -rf /tmp/ipi_h2o.E0_$E0.traj_$traj
rm -rf /tmp/ipi_h2o-dipole_$E0.traj_$traj
#Run noneqm-traj script. Similar to i-pi but requires two inputs: i-pi input and spectra-related input.
./noneqm-traj.py input.xml -e 0.1 -t 40 &> output &
sleep 5
lmp < in.lmp &> out_lammps &
i-pi-driver -u -h h2o-dipole_$E0.traj_$traj -m water_dip_pol -o 0 &> out
wait
python ./noneqm-response.py input.xml -e 0.1 -t 40 -b dip,pol -p 0,0