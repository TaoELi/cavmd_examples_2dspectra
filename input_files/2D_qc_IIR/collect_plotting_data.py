import numpy as np
import glob
import os

output_folder = "./final_data/"
directory = os.path.dirname(output_folder)
if directory and not os.path.exists(directory): os.makedirs(directory, exist_ok=True)

data_folder="64_did"
directory1 = os.path.dirname(output_folder+data_folder)
if directory1 and not os.path.exists(directory1): os.makedirs(directory1, exist_ok=True)

patterns = ["2d-qc-raman.out"]

for e0 in [7]:
    directory2 = os.path.dirname(output_folder+data_folder+f"/{e0}e-4/")
    if directory2 and not os.path.exists(directory2): os.makedirs(directory2, exist_ok=True)
    for name in patterns:
        for i in range(1,11):
            if i == 1 : data = np.loadtxt("./"+data_folder+f"/{i}_{e0}e-4/result/final_results/"+name)
            else : data += np.loadtxt("./"+data_folder+f"/{i}_{e0}e-4/result/final_results/"+name)
        data /= 10
        np.savetxt(output_folder+data_folder+f"/{e0}e-4/"+name, data)
