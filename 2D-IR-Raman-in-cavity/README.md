2D IR-Raman spectrum of liquid water under VSC 
====================================

Author:
`Xinwei Ji <astolfo@udel.edu>`

This example demonstrates the use of equilibrium-nonequilibrium CavMD within the framework of i-PI to simulate 2D IR-Raman spectrum of liquid water under vibrational strong coupling (VSC). 
For theory, see

1. [TBD](TBD)

Also see the [2D IR-Raman spectrum without VSC](https://github.com/i-pi/i-pi/tree/main/demos/2D-IR-Raman).

The code provided here can be used to reproduce the results in reference.

Briefly, the goal is to compute the response function R(t_1, t_2) according to fig 1 and eq 13 of reference, for which we need a backward equilibrium trajectory and two forward nonequilibrium (kicked) trajectories starting from an initial condition sampled according to the thermal distribution. In practice, the equilibrium part is computed by running a single long trajectory and nonequilibrium trajectories are started from points along the equilibrium trajectory. We then average over all "nonequilibrium events". The nonequilibrium part of the calculation is implemented in the enclosed script noneqm-traj.py and the computation of the response function is implemented in noneqm-response.py. 

This example runs on a periodic box of 64 water molecules at 280 K.

This example does not produce useful data, but additional instructions are added at the end of this file for those who wish to run a converged simulation similar to that of reference.

Required software
-----------------

This simulation requires the use of [LAMMPS](https://www.lammps.org/), including the 
`fix_ipi` module and all the modules needed to use a TIP4P water potential.

For the calculation of dipoles, polarizabilities, and dipole derivatives, we use the
water_dip_pol option from the i-PI f90 driver.


Usage
-----

Run all trajectories and compute the response function using the provided run.sh script 

```bash
./run.sh
```

This takes about 2 minutes to finish.

The output is found in dip_pol_neq_2d.dat file. It contains R(t_1, t_2) for 10 fs along both axes.

Step-by-step explanation
-------------------------

**Equilibrium part:**

The equilibrium trajectory that we run is a MD trajectory and also computes properties (dipole moment, polarizability, and dipole gradient) by using a driver that sends these to i-PI inside the extras string of the socket (option 'water_dip_pol' of i-PI f90 driver). To print these out, we added
```xml
<trajectory filename='dip' stride='4' extra_type='dipole'> extras </trajectory>
<trajectory filename='pol' stride='4' extra_type='polarizability'> extras </trajectory>
<trajectory filename='dip_der' stride='20' extra_type='dipole_derivative'> extras </trajectory>
```
The strides of dipole and polarizability defines the timestep for our t_1 and t_2 times in the response function. Stride of the dipole derivative defines the steps at which we will initiate nonequilibrium trajectories (see later). To run nonequilibrium trajectories, we also need their starting positions and momenta, so we also need to print out checkpoint files every 20 steps:
```xml
<checkpoint stride='20' filename='chk' overwrite='false'/>
```
To tell i-PI to use both lammps (for forces come from both molecules and photons) and f90 driver (for properties), we need to add both of these as ffcavphsocket and "forces" (note: water_dip_pol model returns zero force, so it does not affect dynamics). In run.sh, the equilibrium trajectory is then run using:
```bash
i-pi input.xml &> output &
sleep 10
lmp < in.lmp &> out_lammps &
i-pi-driver -u -h h2o-dipole -m water_dip_pol -o 1 &> out
wait
```
Note that i-pi-driver option water_dip_pol accepts one argument that is 1 if we want to compute the gradient of the dipole moment or 0 if not.

To setup the CavMD simulation, we need to define ffcavphsocket part in xml file for lammps.
```xml
<ffcavphsocket name='lammps' mode='unix' pbc='True'>
<address> h2o.E0_7e-4.traj_1 </address>
<latency> 0.01 </latency>
<apply_photon> True </apply_photon>
<E0> 7e-4 </E0>
<omega_c units='inversecm'> 3550.0 </omega_c>
<charge_array> [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] </charge_array>
</ffcavphsocket>
```
and dipole f90 dirver
```xml
<ffcavphsocket name='dipole' mode='unix' pbc="false">
<address> h2o-dipole_7e-4.traj_1 </address>
<latency> 0.01 </latency>
<apply_photon> True </apply_photon>
<E0> 7e-4 </E0>
<omega_c units='inversecm'> 3550.0 </omega_c>
<charge_array> [] </charge_array>
</ffcavphsocket>
```
where apply_photon control whether the system is coupled with photons at freqency omega_c, and E0 denotes the effective coupling strength between light and matter. The charge_array is a array defining the charge on each atom envolved in CavMD propagation. 

If we want employ the DID dipole and dipole derivatives calculated from f90 dipole driver in the CavMD propagation, we need to set a all-zero charge_array which tells the ffcavphsocket do not double count the dipole and dipole derivatives calculated in lammps part, and a empty charge_array [] in dipole part.

In another scenario, if we want employ the point charge dipole and dipole derivatives in the CavMD propagation, we need to set a non-zero charge_array in lammps part, and a all-zero charge_array in dipole part.


**Nonequilibrium part**:

For every checkpoint file we created in the previous step, we will restart two trajectories: One trajectory will be started with (q, p - epsilon * dip_der / 2) 
and another trajectory will use (q, p + epsilon * dip_der / 2), where q and p are the position and momentum vectors.
Here, epsilon is a parameter that defines the magnitude of the momentum kick (simulating an instantaneous interaction with an external filed, 
so epsilon is proportional to the field's amplitude). dip_der is the coordinate gradient of the x-component of the dipole moment, which is read from the outputs of the
equilibrium trajectory calculation.

To reduce the amount of time spent on initializing i-PI, sockets, drivers, we use a script that performs above calculation for all trajectories
without closing the sockets. This script can be run in a similar way as the main i-PI executable:
```bash
./noneqm-traj.py input.xml -e 0.1 -t 40 &> output &
sleep 10
$lmp < in.lmp &> out_lammps &
i-pi-driver -u -h h2o-dipole -m water_dip_pol -o 0 &> out
wait
```
The noneqm-traj.py script takes as input the original input.xml file (which is used to get information like simulation name, timestep, temperature, integrators, ...)
and two paramaters:
1) epsilon (option -e, here set to 0.1)
2) number of time steps for nonequilibrium dynamics (option -t, here set to 40).
Note that the number of time steps provided by the -t option will also determine the maximum t_2 time for R(t_1, t_2).
Here, we can turn off the computation of dipole gradients ('-o 0' option for the i-pi-driver).

To compute the response function later, we will need backward dynamics along t_1, which means that we have no use of running nonequilibrium trajectories from 
the first few checkpoint files. The script is set to ensure that we have same maximum t_1 and t_2 times (in our case 40 steps or 40 x 0.25 fs = 10 fs). This
corresponds to the checkpoint file simulation.chk_2 (printed at step 40 of equilibrium dynamics). In total, we will run 9 chk files x 2 trajectories = 18 nonequilibrium trajectories.
The properties of these trajectories are printed in simulation_noneqm.dip/pol_X (X is the bead number).

noneqm-traj.py script also accepts optional parameters -i (--tfirst) and -f (--tlast) if the user wants to specify the starting and ending (included) checkpoint files (see ./noneqm-traj.py --help). This might be useful if the simulation cost requires us to split the simulation into chunks.

As a sanity check, you can see that the dipole moment (e.g., of bead 0) is same at step 40 of equilibrium trajectory (line 22 of simulation.dip_0), at step 0 of the first nonequilibrium trajectory (line 2 of simulation_noneqm.dip_0), and at step 0 of the second nonequilibrium trajectory (line 24 of simulation_noneqm.dip_0). This is because at time zero the nonequilibrium trajectories have the same position. Because of the kicks, all these trajectories (and dipoles/polarizabilities) become different after the first step.

**Response function**:

noneqm-response.py uses the outputs of above calculations to evaluate R(t_1, t_2). In this example, the second interaction is with a dipole moment (this was determined by using a dipole gradient to perform the momentum kicks). We can still choose the type of spectroscopy we want to simulate by setting the operators that we wish to correlate. 

To see the usage of noneqm-response.py, run
```bash
./noneqm-response.py --help
```
First two arguments -e and -t should coincide with those used for noneqm-traj.py above. -t sets the maximum t_1 and t_2 times. The step is determined by the output stride of relevant quantities (dip, pol) in the dynamics above. For us, this was 4 dynamics time steps or 1 fs. With -t 40, we have 10 fs of response function along t_1 and 10 fs along t_2.
If noneqm-traj.py was run using --tfirst and/or --tlast, same options should be used here.

The main additional options here are -b and -p:
--op or -b option can be set to 'dip, pol', which will compute the IR-IR-Raman (or IIR) response function. This is default. 
Since the dipole moment operators are vectors and polarizability a tensor, the response function is in general a rank-4 tensor [Eq.13 in ref]. The photon couples with the molecules along x-axis, and the dipole gradient is already computed along x-axis, so we know we will have R_{axcd}, where a=x,y,z will correspond to the dipole component along negative t_1 axis and cd = xx, xy, ... correspond to the component of the polarizability along the nonequilibrium trajectories (time axis t_2). This can be set by -p or --field_pol option, whose default is '0, 0' (this corresponds to the xxxx component of the response function).

We can now compute the response function as shown under "Usage", which is equivalent to
```bash
python ./noneqm-response.py input.xml -e 0.1 -t 40 -b dip,pol -p 0,0
```

How to obtain meaningful results?
---------------------------------

1) Our example here uses only NVE trajectories. To sample from the thermal distribution, we need to perform this type of calculation starting from many initial conditions sampled from the appropriate ensemble (e.g., NVT).
2) Here we run a very short equilibrium trajectory, which allows us to average over only 9 nonequilibrium events. Also, taking only 20 steps as the stride for nonequilibrium trajectory calculations makes little sense. A more reasonable value for liquid water would be around 1000 steps.
3) To compute a meaningful 2D spectra of liquid water, t_1 and t_2 times should be on the order of 100 fs - few ps.

For example, ref used:
1) 5120 initial samples from an equilibrated NVT calculation as starting points for the equilibrium-nonequilibrium trajectory runs.
2) Equilibrium trajectories were 25 ps (100,000 steps) long and checkpoint/dipole-derivative stride was 1000 steps (250 fs), resulting in about 100 nonequilibrium events.
3) Nonequilibrium trajectories were 1000 steps (250 fs) long.

An example of a reasonably converged calculation for 250 fs long response function is available in file dip_pol_neq_2d_in_cavity.converged. Corresponding input.xml and run.sh files can be found in inputs_for_converged_calculation directory.


How to compute a spectrum from the response function?
-----------------------------------------------------

See spectrum.ipynb Jupyter notebook, which uses a 2D sine transform to compute the spectrum from the response function in dip_pol_neq_2d.converged.
