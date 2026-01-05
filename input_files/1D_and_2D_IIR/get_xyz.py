# This script is used to capture many important data from each xyz trajectory

import numpy as np
import time
import MDAnalysis as mda    

class MD_Analysis:
    def __init__(self, xyz_filename, dtfs = 2, nframe_max=114514):
        self.xyz_filename = xyz_filename
        self.dtfs = dtfs
        self.dtau = dtfs * 1e-15 / 2.418884326e-17
        self.nframe_max = nframe_max
        self.natoms = 0
        self.labels = []
        self.traj = []
        self.nmolecule = 0
        self.load_xyz(self.xyz_filename)
        # After read xyz file, now we calculate different properties
        
    def load_xyz(self, xyz_filename):
        time_start = time.time()
        with open(xyz_filename, 'r') as myfile:
            self.natoms = int(myfile.readline().strip())
        data_raw = mda.Universe(xyz_filename)
        traj = data_raw.trajectory
        nframes = len(traj)
        natoms = np.shape(np.array(traj[0]))[0]
        trajectory_data = np.zeros((natoms, 3, nframes))
        for ts in range(nframes):
            #trajectory_data[:,:,ts] = np.array(traj[ts])
            trajectory_data[:,:,ts] = traj[ts]._pos.copy()
        print("There are %d frames in %s" %(nframes, xyz_filename))
        self.traj = trajectory_data[:self.natoms,:,:]
        print(np.shape(self.traj))
        atom = data_raw.atoms
        text = f'{self.natoms}\n'
        for tstep in range(1,513):
            with open(f'./init_water_xyz_nvt/water_64.xyz_{tstep}', 'w') as output:
                print(text, file=output)
                for i in range(len(atom)): print(f'       {atom[i].element} {trajectory_data[i,0,tstep]:12.5e} {trajectory_data[i,1,tstep]:12.5e} {trajectory_data[i,2,tstep]:12.5e}', file=output)
        time_end = time.time()
        print(f'read traj cost time = {time_end-time_start:.2f} s')

if __name__ == "__main__":
    a = MD_Analysis(xyz_filename="./simu.pos_0.xyz")
