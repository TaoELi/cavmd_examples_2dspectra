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

patterns = ["simu_*.txt", "out_*_neq_2d.dat_*"]
data5 = obtain_avg_data(path=f"{input}/noneq-qc-raman", pattern=patterns[1])
np.savetxt(f'{output}/2d-qc-raman.out', data5)
