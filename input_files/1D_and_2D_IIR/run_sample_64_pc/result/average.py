import numpy as np
import glob

input=None
output=None

def obtain_avg_data(path, pattern="simu_*.vac.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories" %(len(filenames)))
    data = np.loadtxt(filenames[0])
    for filename in filenames[1:]:
        #print(filename)
        data += np.loadtxt(filename)
    data /= float(len(filenames))
    return data

patterns = ["simu_*.txt", "dip_*_neq_2d.dat_*"]
data0 = obtain_avg_data(path=f"{input}/phdac_simple_charge", pattern=patterns[0])
np.savetxt(f'{output}/phdac.out', data0)
data1 = obtain_avg_data(path=f"{input}/dac_simple_charge", pattern=patterns[0])
np.savetxt(f'{output}/dacpc.out', data1)
data2 = obtain_avg_data(path=f"{input}/dac_dipole_surface", pattern=patterns[0])
np.savetxt(f'{output}/dacdid.out', data2)
data3 = obtain_avg_data(path=f"{input}/pac_iso_dipole_surface", pattern=patterns[0])
np.savetxt(f'{output}/pacdid_iso.out', data3)
data4 = obtain_avg_data(path=f"{input}/pac_aniso_dipole_surface", pattern=patterns[0])
np.savetxt(f'{output}/pacdid_aniso.out', data4)
data5 = obtain_avg_data(path=f"{input}/noneq-ir-raman", pattern=patterns[1])
np.savetxt(f'{output}/2d-ir-raman.out', data5)
