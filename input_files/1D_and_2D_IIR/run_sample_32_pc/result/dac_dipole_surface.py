# This script is used to capture many important data from each xyz trajectory

import numpy as np
from scipy import signal
import os, sys
from itertools import islice
from scipy import fftpack
import glob
import time
import MDAnalysis as mda    

def smooth(x,window_len=11,window='hamming'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if window_len<3:
        return x

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[window_len//2-1:-window_len//2]

class MD_Analysis:
    def __init__(self, xyz_filename, dtfs = 1, nframe_max=114514):
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
        data = np.loadtxt(xyz_filename)
        self.dipole_water_traj = np.transpose(data)

    def cacl_dipole_traj(self):
        print("Calculating the autocorrelation function...")
        self.dacf_x = self.auto_correlation_function_simple(self.dipole_water_traj[0,:])
        self.dacf_y = self.auto_correlation_function_simple(self.dipole_water_traj[1,:])
        self.dacf_z = self.auto_correlation_function_simple(self.dipole_water_traj[2,:])
        self.dacf_tot = self.dacf_x + self.dacf_y + self.dacf_z
        self.dacf_time_fs = np.linspace(0.0, self.dtfs*(self.dacf_x.size -1), self.dacf_x.size)
        print("Calculating the FFT of autocorrelation function")
        self.dacf_x_freq, self.dacf_x_sp = self.fft3(self.dacf_x)
        self.dacf_y_freq, self.dacf_y_sp = self.fft3(self.dacf_y)
        self.dacf_z_freq, self.dacf_z_sp = self.fft3(self.dacf_z)
        self.dacf_tot_freq, self.dacf_tot_sp = self.fft3(self.dacf_tot)

    def cacl_pol_traj(self):
        print("Calculating the autocorrelation function...")
        self.pacf_x = self.auto_correlation_function_simple(self.dipole_water_traj[0,:])
        self.pacf_y = self.auto_correlation_function_simple(self.dipole_water_traj[3,:])
        self.pacf_z = self.auto_correlation_function_simple(self.dipole_water_traj[6,:])
        self.pacf_tot = self.pacf_x + self.pacf_y + self.pacf_z
        self.pacf_time_fs = np.linspace(0.0, self.dtfs*(self.pacf_x.size -1), self.pacf_x.size)
        print("Calculating the FFT of autocorrelation function")
        self.pacf_x_freq, self.pacf_x_sp = self.fft3(self.pacf_x)
        self.pacf_y_freq, self.pacf_y_sp = self.fft3(self.pacf_y)
        self.pacf_z_freq, self.pacf_z_sp = self.fft3(self.pacf_z)
        self.pacf_tot_freq, self.pacf_tot_sp = self.fft3(self.pacf_tot)

    def auto_correlation_function_fft(self, x):
        corr = signal.fftconvolve(x, x[::-1], mode='same')
        corr = corr[corr.size // 2: ]
        return corr / corr[0] * np.mean(x * x)

    def auto_correlation_function_simple(self, x):
        n = x.size
        if n % 2 == 0:
            x_shifted = np.zeros(n*2)
        else:
            x_shifted = np.zeros(n*2-1)
        x_shifted[n//2 : n//2+n] = x
        # Convolute the shifted array with the flipped array, which is equivalent to performing a correlation
        autocorr_full = (signal.fftconvolve(x_shifted, x[::-1], mode='same')[-n:]/ np.arange(n, 0, -1))
        # Truncate the autocorrelation array
        autocorr = autocorr_full[0:n//2]
        return autocorr

    def fft(self, x):
        #sp = np.fft.fft(x)
        #freq_au = 2.0 * np.pi * np.fft.fftfreq(np.size(x), self.dtau)
        # Because dt has the unit of fs, I need to transform fs^{-1} to cm^{-1}
        #freq_cminverse = freq_au * 219474.63
        #return freq_cminverse[0:sp.size//2], sp[0:sp.size//2]
        lineshape = fftpack.dct(x, type=1)
        freq_au = np.linspace(0, 0.5/self.dtfs * 1e15, len(x))
        # Because dt has the unit of fs, I need to transform fs^{-1} to cm^{-1}
        freq_cminverse = freq_au / (100.0 * 299792458.0)
        return freq_cminverse, lineshape

    def fft3(self, x ):
        # Adding zeros to the end of x
        lineshape = fftpack.dct(x, type=1)
        freq_au = np.linspace(0, 0.5/self.dtfs * 1e15, len(x))
        # Because dt has the unit of fs, I need to transform fs^{-1} to cm^{-1}
        freq_cminverse = freq_au / (100.0 * 299792458.0)
        # Calculate spectra
        #field_description =  freq_au**2
        field_description =  freq_au**2
        spectra = lineshape * field_description
        return freq_cminverse, spectra
        #return freq_cminverse[0:spectra.size//2], spectra[0:spectra.size//2]

    def output_dipole_autocorrelation(self):
        local_filename = "%s.dac.txt" %self.xyz_filename
        if os.path.isfile(local_filename + "dd"):
            print("Have calculated dipole autocorrelation for %s, skipping..." %self.xyz_filename)
        else:
            self.cacl_dipole_traj()
            # output data
            data = np.zeros((np.size(self.dacf_x), 10))
            data[:, 0] = self.dacf_time_fs
            data[:, 1] = self.dacf_x
            data[:, 2] = self.dacf_y
            data[:, 3] = self.dacf_z
            data[:, 4] = self.dacf_tot
            data[:, 5] = self.dacf_x_freq
            data[:, 6] = smooth(self.dacf_x_sp)
            data[:, 7] = smooth(self.dacf_y_sp)
            data[:, 8] = smooth(self.dacf_z_sp)
            data[:, 9] = smooth(self.dacf_tot_sp)
            comments = "# dacf_time_fs, dacf_x, dacf_y, dacf_z, dacf_tot, freq, sp_x, sp_y, sp_z, sp_tot"
            np.savetxt(local_filename, data, comments=comments)

    def output_pol_autocorrelation(self):
        local_filename = "%s.pac.txt" %self.xyz_filename
        if os.path.isfile(local_filename + "dd"):
            print("Have calculated dipole autocorrelation for %s, skipping..." %self.xyz_filename)
        else:
            self.cacl_pol_traj()
            # output data
            data = np.zeros((np.size(self.pacf_x), 10))
            data[:, 0] = self.pacf_time_fs
            data[:, 1] = self.pacf_x
            data[:, 2] = self.pacf_y
            data[:, 3] = self.pacf_z
            data[:, 4] = self.pacf_tot
            data[:, 5] = self.pacf_x_freq
            data[:, 6] = smooth(self.pacf_x_sp)
            data[:, 7] = smooth(self.pacf_y_sp)
            data[:, 8] = smooth(self.pacf_z_sp)
            data[:, 9] = smooth(self.pacf_tot_sp)
            comments = "# dacf_time_fs, dacf_x, dacf_y, dacf_z, dacf_tot, freq, sp_x, sp_y, sp_z, sp_tot"
            np.savetxt(local_filename, data, comments=comments)

if __name__ == "__main__":
    paths = sys.argv[1:]
    for i, path in enumerate(paths):
        filenames = glob.glob("%s/simu_*.dip_0" %path)
        print(len(filenames))
        for filename in filenames:
            a = MD_Analysis(xyz_filename=filename)
            a.output_dipole_autocorrelation()
            #a.output_pol_autocorrelation()
