import matplotlib.pyplot as plt
import numpy as np
import columnplots as clp
import scipy.fft as fftpack
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import matplotlib as mpl

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

def plot_qc():
    axes = clp.initialize(2, 1, width=4.3, height=4.3*0.618*1.1*2, LaTeX=True, sharex=True, sharey=True, fontsize=12)
    e0list = [4,7]
    ys_sc_sc, ys_ds_sc = [], []
    for e0 in range(len(e0list)):
        data = np.loadtxt(f'./final_results/qc/qc_pc_{e0list[e0]}e-4.out')
        x = data[:,4]
        interval = np.where(x<4500)[0]
        y = data[:,7]/2e28
        ymax = abs(y[interval]).max()
        y = y / ymax * 0.8 + e0
        ys_sc_sc.append(y[interval])
        data = np.loadtxt(f'./final_results/qc/qc_did_{e0list[e0]}e-4.out')
        y = data[:,7]/2e28
        ymax = abs(y[interval]).max()
        y = y / ymax * 0.8 + e0
        ys_ds_sc.append(y[interval])
    
    xs = [x[interval]] * 2
    clp.plotone(xs, ys_sc_sc, axes[0], lw=1, 
                xlim=[0, 4500], ylim=[-0.05,2], 
                xlabel=None, ylabel="IR intensity [arb. units]", 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    clp.plotone(xs, ys_ds_sc, axes[1], lw=1, 
                xlim=[0, 4500], ylim=[-0.05,2], 
                xlabel=r"frequency [cm$^{-1}$]", ylabel="IR intensity [arb. units]", 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)

    axes[0].axhline(y=0, linestyle='dotted', alpha=0.3, linewidth=1)
    axes[0].axhline(y=1, linestyle='dotted', alpha=0.3, linewidth=1)
    axes[1].axhline(y=0, linestyle='dotted', alpha=0.3, linewidth=1)
    axes[1].axhline(y=1, linestyle='dotted', alpha=0.3, linewidth=1)

    cadjust_colors = [plt.cm.hot(i) for i in np.linspace(0, 0.6, 4)]
    name_list = ["(a) PC in CavMD", "(b) DID in CavMD"]

    for j in range(2):
        for e0 in range(2):
            if e0 == 0 : axes[j].text(1000, e0+0.3, r"$I_{\rm{c}}^{\rm{IR}}, \ \widetilde{\varepsilon}=4\times 10^{-4}$ a.u.", fontsize=10, color=cadjust_colors[e0])
            elif e0 == 1 : axes[j].text(1000, e0+0.3, r"$I_{\rm{c}}^{\rm{IR}}, \ \widetilde{\varepsilon}=7\times 10^{-4}$ a.u.", fontsize=10, color=cadjust_colors[e0])
        axes[j].text(100, 1.8, name_list[j], fontsize=10, color="k")
        axes[j].set_yticks([0,1,2])

    plt.rcParams["axes.axisbelow"] = False
    clp.adjust(savefile=f"./fig/1d-qc.png")

def plot_IR():
    axes = clp.initialize(2, 2, width=4.3*2, height=4.3*0.618*1.1*2, LaTeX=True, sharex=True, sharey=True, fontsize=12)
    e0list = [4,7]
    data_outside_sc = np.loadtxt(f'./final_results/ir/dac_did_pc_0e-4.out')
    data_outside_ds = np.loadtxt(f'./final_results/ir/dac_did_did_0e-4.out')
    x, y_outside_sc = data_outside_sc[:,5], (data_outside_sc[:,6] + data_outside_sc[:,7])/2e28
    x, y_outside_ds = data_outside_ds[:,5], (data_outside_ds[:,6] + data_outside_ds[:,7])/2e28
    interval = np.where(x<4500)[0]
    norm_interval = np.where(x<1000)[0]
    norm_factor = 0.07
    xs = [x[interval]] * 3
    y_outside_sc_max = abs(y_outside_sc[norm_interval]).max()
    y_outside_ds_max = abs(y_outside_ds[norm_interval]).max()
    ys_sc_sc, ys_ds_sc = [y_outside_sc[interval] / y_outside_sc_max * norm_factor], [y_outside_sc[interval] / y_outside_sc_max * norm_factor]
    ys_sc_ds, ys_ds_ds = [y_outside_ds[interval] / y_outside_ds_max * norm_factor], [y_outside_ds[interval] / y_outside_ds_max * norm_factor]

    for e0 in range(len(e0list)):
        data = np.loadtxt(f'./final_results/ir/dac_pc_pc_{e0list[e0]}e-4.out')
        y = (data[:,6] + data[:,7])/2e28
        ymax = abs(y[norm_interval]).max()
        y = y / ymax * norm_factor + e0 + 0.5
        ys_sc_sc.append(y[interval])
        data = np.loadtxt(f'./final_results/ir/dac_pc_did_{e0list[e0]}e-4.out')
        y = (data[:,6] + data[:,7])/2e28
        ymax = abs(y[norm_interval]).max()
        y = y / ymax * norm_factor + e0 + 0.5
        ys_sc_ds.append(y[interval])
        data = np.loadtxt(f'./final_results/ir/dac_did_pc_{e0list[e0]}e-4.out')
        y = (data[:,6] + data[:,7])/2e28
        ymax = abs(y[norm_interval]).max()
        y = y / ymax * norm_factor + e0 + 0.5
        ys_ds_sc.append(y[interval])
        data = np.loadtxt(f'./final_results/ir/dac_did_did_{e0list[e0]}e-4.out')
        y = (data[:,6] + data[:,7])/2e28
        ymax = abs(y[norm_interval]).max()
        y = y / ymax * norm_factor + e0 + 0.5
        ys_ds_ds.append(y[interval])

    clp.plotone(xs, ys_sc_sc, axes[0,0], lw=1, 
                xlim=[0, 4500], ylim=[-0.2,2.5], 
                xlabel=None, ylabel="IR intensity [arb. units]", 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    clp.plotone(xs, ys_sc_ds, axes[1,1], lw=1, 
                xlim=[0, 4500], ylim=[-0.2,2.5], 
                xlabel=r"frequency [cm$^{-1}$]", ylabel=None, 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    clp.plotone(xs, ys_ds_sc, axes[1,0], lw=1, 
                xlim=[0, 4500], ylim=[-0.2,2.5], 
                xlabel=r"frequency [cm$^{-1}$]", ylabel="IR intensity [arb. units]", 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    clp.plotone(xs, ys_ds_ds, axes[0,1], lw=1, 
                xlim=[0, 4500], ylim=[-0.2,2.5], 
                xlabel=None, ylabel=None, 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    axes[1,1].annotate('', xy=(3590, 0.55), xytext=(3590, 0.8), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[1,1].annotate('', xy=(3590, 1.55), xytext=(3590, 1.8), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[1,0].annotate('', xy=(3600, 0.55), xytext=(3600, 0.8), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[1,0].annotate('', xy=(3600, 1.55), xytext=(3600, 1.8), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    
    cadjust_colors = [plt.cm.hot(i) for i in np.linspace(0, 0.6, 3)]
    name_list = ["(a) PC-PC", "(b) DID-DID", "(c) DID-PC", "(d) PC-DID"]
    for j in range(4):
        xidx, yidx = j // 2, j % 2
        for e0 in range(3):
            if e0 == 0 : axes[xidx, yidx].text(1500, 0.2, "outside cavity", fontsize=10, color=cadjust_colors[e0])
            elif e0 == 1 : axes[xidx, yidx].text(1500, e0-0.2, r"$\widetilde{\varepsilon}=4\times 10^{-4}$ a.u.", fontsize=10, color=cadjust_colors[e0])
            else: axes[xidx, yidx].text(1500, e0-0.2, r"$\widetilde{\varepsilon}=7\times 10^{-4}$ a.u.", fontsize=10, color=cadjust_colors[e0])
        axes[xidx, yidx].text(200, 2.2, name_list[j], fontsize=10, color="k")

    plt.rcParams["axes.axisbelow"] = False
    clp.adjust(savefile=f"./fig/1d-ir.png")

def plot_IR_size():
    axes = clp.initialize(2, 2, width=4.3*2, height=4.3*0.618*1.1*2, LaTeX=True, sharex=True, sharey=True, fontsize=12)
    nlist = [32,64,128]
    ys_sc_sc, ys_ds_sc, ys_sc_ds, ys_ds_ds = [], [], [], []

    for n in nlist:
        data = np.loadtxt(f'./final_results/ir_size/dac_7e-4_{n}_pc_pc.out')
        x, y = data[:,5], (data[:,6] + data[:,7])/2e28
        interval = np.where(x<4500)[0]
        xs = [x[interval]] * 3
        ymax = abs(y[np.where(x<1500)[0]]).max()
        y = y / ymax
        ys_sc_sc.append(y[interval])
        data = np.loadtxt(f'./final_results/ir_size/dac_7e-4_{n}_pc_did.out')
        y = (data[:,6] + data[:,7])/2e28
        ymax = abs(y[np.where(x<1500)[0]]).max()
        y = y / ymax
        ys_sc_ds.append(y[interval])
        data = np.loadtxt(f'./final_results/ir_size/dac_7e-4_{n}_did_pc.out')
        y = (data[:,6] + data[:,7])/2e28
        ymax = abs(y[np.where(x<1500)[0]]).max()
        y = y / ymax
        ys_ds_sc.append(y[interval])
        data = np.loadtxt(f'./final_results/ir_size/dac_7e-4_{n}_did_did.out')
        y = (data[:,6] + data[:,7])/2e28
        ymax = abs(y[np.where(x<1500)[0]]).max()
        y = y / ymax
        ys_ds_ds.append(y[interval])

    labels = [r"$N_{\rm{simu}}=%d$"%ni for ni in nlist]
    clp.plotone(xs, ys_sc_sc, axes[0,0], lw=1, 
                xlim=[0, 4500], ylim=[-0.2,20], 
                labels=labels,
                xlabel=None, ylabel="IR intensity [arb. units]", 
                showlegend=True, legendloc=(0.3,0.4), legendFontSize=9,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    clp.plotone(xs, ys_sc_ds, axes[1,1], lw=1, 
                xlim=[0, 4500], ylim=[-0.2,20], 
                xlabel=r"frequency [cm$^{-1}$]", ylabel=None, 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    axins01 = axes[1,1].inset_axes([0.15, 0.3, 0.35, 0.4])  # 位置/大小
    clp.plotone(xs, ys_sc_ds, axins01, lw=1, 
                xlim=[3400, 3800], ylim=[-0.05,0.5], 
                xlabel=None, ylabel=None, 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    axes[1,1].indicate_inset_zoom(axins01, edgecolor='gray')
    clp.plotone(xs, ys_ds_sc, axes[1,0], lw=1, 
                xlim=[0, 4500], ylim=[-0.2,20], 
                xlabel=r"frequency [cm$^{-1}$]", ylabel="IR intensity [arb. units]", 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    axins10 = axes[1,0].inset_axes([0.15, 0.3, 0.35, 0.4])  # 位置/大小
    clp.plotone(xs, ys_ds_sc, axins10, lw=1, 
                xlim=[3400, 3800], ylim=[-0.02,0.2], 
                xlabel=None, ylabel=None, 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    axes[1,0].indicate_inset_zoom(axins10, edgecolor='gray')
    clp.plotone(xs, ys_ds_ds, axes[0,1], lw=1, 
                xlim=[0, 4500], ylim=[-0.2,20], 
                xlabel=None, ylabel=None, 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)

    cadjust_colors = [plt.cm.hot(i) for i in np.linspace(0, 0.6, 11)]
    name_list = ["(a) PC-PC", "(b) DID-DID", "(c) DID-PC", "(d) PC-DID"]
    for j in range(4):
        xidx, yidx = j // 2, j % 2
        axes[xidx, yidx].set_yticks([1,5,10,15,20])
        axes[xidx, yidx].text(100, 18, name_list[j], fontsize=10, color="k")
        axes[xidx, yidx].axvline(x=650, ymin=0, ymax=0.05, color='grey', linestyle='--', alpha=0.5)
        axes[xidx, yidx].axhline(y=1, xmin=0, xmax=0.15, color='grey', linestyle='--', alpha=0.5)

    plt.rcParams["axes.axisbelow"] = False
    clp.adjust(savefile=f"./fig/1d-ir_size_dependence.png")

def plot_Raman():
    axes = clp.initialize(2, 2, width=4.3*2, height=4.3*0.618*1.1*2, LaTeX=True, fontsize=12, sharey=True, sharex=True)
    e0list=[4,7]
    data_outside_iso = np.loadtxt(f'./final_results/raman/did_0e-4/pac_iso.out')
    x, y = data_outside_iso[:,2], data_outside_iso[:,3]/1e28
    interval = np.where((10<x)&(x<4500))[0]
    ymax = abs(y[interval]).max()
    y = y / ymax * 0.8
    xs = [x[interval]] * 3
    ys_sc_ds_iso = [y[interval]]
    ys_ds_ds_iso = [y[interval]]
    data_outside_aniso = np.loadtxt(f'./final_results/raman/did_0e-4/pac_aniso.out')
    x, y = data_outside_aniso[:,2], data_outside_aniso[:,3]/1e28
    ymax = abs(y[interval]).max()
    y = y / ymax * 0.8
    ys_sc_ds_aniso = [y[interval]]
    ys_ds_ds_aniso = [y[interval]]
    
    for e0 in range(len(e0list)):
        datasc_iso = np.loadtxt(f'./final_results/raman/pc_{e0list[e0]}e-4/pac_iso.out')
        ysc = datasc_iso[:,3]/1e28
        ymax = abs(ysc[interval]).max()
        ysc = ysc / ymax * 0.8 + e0 + 1
        ys_sc_ds_iso.append(ysc[interval])
        datads_iso = np.loadtxt(f'./final_results/raman/did_{e0list[e0]}e-4/pac_iso.out')
        yds = datads_iso[:,3]/1e28
        ymax = abs(yds[interval]).max()
        yds = yds / ymax * 0.8 + e0 + 1
        ys_ds_ds_iso.append(yds[interval])
        datasc_aniso = np.loadtxt(f'./final_results/raman/pc_{e0list[e0]}e-4/pac_aniso.out')
        ysc = datasc_aniso[:,3]/1e28
        ymax = abs(ysc[interval]).max()
        ysc = ysc / ymax * 0.8 + e0 + 1
        ys_sc_ds_aniso.append(ysc[interval])
        datads_aniso = np.loadtxt(f'./final_results/raman/pc_{e0list[e0]}e-4/pac_aniso.out')
        yds = datads_aniso[:,3]/1e28
        ymax = abs(yds[interval]).max()
        yds = yds / ymax * 0.8 + e0 + 1
        ys_ds_ds_aniso.append(yds[interval])

    clp.plotone(xs, ys_sc_ds_iso, axes[0,0], lw=1.5, 
                xlim=[0, 4500], ylim=[-0.2,3], 
                xlabel=None, 
                ylabel="Raman intensity [arb. units]", 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    clp.plotone(xs, ys_ds_ds_iso, axes[0,1], lw=1.5, 
            xlim=[0, 4500], ylim=[-0.2,3], 
            xlabel=None, 
            ylabel=None, 
            showlegend=False,
            colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    clp.plotone(xs, ys_sc_ds_aniso, axes[1,0], lw=1.5, 
                xlim=[0, 4500], ylim=[-0.2,3], 
                xlabel=r"frequency [cm$^{-1}$]", 
                ylabel="Raman intensity [arb. units]", 
                showlegend=False,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    clp.plotone(xs, ys_ds_ds_aniso, axes[1,1], lw=1.5, 
            xlim=[0, 4500], ylim=[-0.2,3], 
            xlabel=r"frequency [cm$^{-1}$]", 
            ylabel=None, 
            showlegend=False,
            colorMap=plt.cm.hot, colorMap_endpoint=0.6)
    
    axes[0,0].annotate('', xy=(3344, 1.25), xytext=(3344, 1.5), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[0,1].annotate('', xy=(3289, 1.25), xytext=(3289, 1.5), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[0,0].annotate('', xy=(3882, 1.25), xytext=(3882, 1.5), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[0,1].annotate('', xy=(3848, 1.25), xytext=(3848, 1.5), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[0,0].annotate('', xy=(3234, 2.15), xytext=(3234, 2.4), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[0,1].annotate('', xy=(3156, 2.15), xytext=(3156, 2.4), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[0,0].annotate('', xy=(4150, 2.15), xytext=(4150, 2.4), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[0,1].annotate('', xy=(4075, 2.15), xytext=(4075, 2.4), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[1,0].annotate('', xy=(3344, 1.25), xytext=(3344, 1.5), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[1,1].annotate('', xy=(3344, 1.25), xytext=(3344, 1.5), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[1,0].annotate('', xy=(3882, 1.25), xytext=(3882, 1.5), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[1,1].annotate('', xy=(3882, 1.25), xytext=(3882, 1.5), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[1,0].annotate('', xy=(3234, 2.15), xytext=(3234, 2.4), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[1,1].annotate('', xy=(3234, 2.15), xytext=(3234, 2.4), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[1,0].annotate('', xy=(4150, 2.15), xytext=(4150, 2.4), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)
    axes[1,1].annotate('', xy=(4150, 2.15), xytext=(4150, 2.4), arrowprops=dict(facecolor='blue', edgecolor='blue', arrowstyle='->', alpha=0.8), fontsize=8)

    cadjust_colors = [plt.cm.hot(i) for i in np.linspace(0, 0.6, 3)]
    for k in range(2):
        for j in range(2):
            for e0 in range(3):
                if e0 == 0 : axes[k,j].text(1200, 0.3, "outside cavity", fontsize=10, color=cadjust_colors[e0])
                elif e0 == 1 : axes[k,j].text(1200, e0+0.3, r"$\widetilde{\varepsilon}=4\times 10^{-4}$ a.u.", fontsize=10, color=cadjust_colors[e0])
                else: axes[k,j].text(1200, e0+0.3, r"$\widetilde{\varepsilon}=7\times 10^{-4}$ a.u.", fontsize=10, color=cadjust_colors[e0])

    axes[0,0].text(100, 2.7, r"(a) $\langle\boldsymbol{\alpha}(0)\boldsymbol{\alpha}(t)\rangle$, PC-DID", fontsize=10, color="k")
    axes[0,1].text(100, 2.7, r"(b) $\langle\boldsymbol{\alpha}(0)\boldsymbol{\alpha}(t)\rangle$, DID-DID", fontsize=10, color="k")
    axes[1,0].text(100, 2.7, r"(c) $\langle\mathrm{Tr}[\boldsymbol{\beta}(0)\boldsymbol{\beta}(t)]\rangle$, PC-DID", fontsize=10, color="k")
    axes[1,1].text(100, 2.7, r"(d) $\langle\mathrm{Tr}[\boldsymbol{\beta}(0)\boldsymbol{\beta}(t)]\rangle$, DID-DID", fontsize=10, color="k")

    plt.rcParams["axes.axisbelow"] = False
    clp.adjust(savefile=f"./fig/1d-raman.png")

def plot_2D_IR_Raman():
    
    fig, axes = clp.initialize(5, 3, width=4.3*3, height=4.3*0.85*5, LaTeX=True, fontsize=12, sharey=True, sharex=True, return_fig_args=True)
    # Useful constants.
    fstoau = 41.34137  # Femtosecond to a.u.
    autocm1 = 219474.63  # a.u. (Hartree energy unit) to wavenumber
    invfstocm1 = autocm1 / fstoau  # fs^-1 to cm^-1 conversion
    # Processing parameters
    dt = 1 * fstoau  # Timestep 1 fs
    nsteps = 250  # Number of steps in Rt
    npad = 1125*2  # Padding for FFT
    nplot = 175*2  # Number of frequency steps to plot.
    dw = 2 * np.pi / (nsteps + npad + 1) / dt * autocm1  # Frequency step for fft
    wmin = 0  # Minimum frequency in the plot
    wmax = nplot * dw  # Maximum frequency in the plot
    # Damping to avoid a hard cut off at long t_1/t_2.
    tau = 1e27
    damping = np.array(
        [
            [np.exp(-(((i + j) * dt / fstoau) ** 12) / tau) for i in range(nsteps + 1)]
            for j in range(nsteps + 1)
        ]
    )
    nlist = [32,64,128] 
    labels = [r"$N_{\rm{simu}}=%d$"%ni for ni in nlist]
    epsilon_list = ["0e-4","4e-4_did","4e-4_pc","7e-4_did","7e-4_pc"]
    text_list = ["outside cavity",r"$\widetilde{\varepsilon}\cdot\sqrt{N_{\rm{simu}}/64}=4 \times 10^{-4}$ a.u.",r"$\widetilde{\varepsilon}\cdot\sqrt{N_{\rm{simu}}/64}=4 \times 10^{-4}$ a.u.",r"$\widetilde{\varepsilon}\cdot\sqrt{N_{\rm{simu}}/64}=7 \times 10^{-4}$ a.u.",r"$\widetilde{\varepsilon}\cdot\sqrt{N_{\rm{simu}}/64}=7 \times 10^{-4}$ a.u."]
    idx = {0:["(a)","(b)","(c)"], 1:["(d)","(e)","(f)"], 2:["(g)","(h)","(i)"], 3:["(j)","(k)","(l)"], 4:["(m)","(n)","(o)"]}
    dipole_list = ["", ", DID-DID", ", PC-DID", ", DID-DID", ", PC-DID"]
    
    for i in range(3):
        for j in range(5):
            # Import Rt.
            Rt = np.loadtxt(f"./final_results/2d_iir_size/2d-ir-raman_{nlist[i]}_{epsilon_list[j]}.out")[:nsteps+1,:nsteps+1]
            # Rt = np.gradient(RtRaw,dt,axis=0,edge_order=2)
            # Pad with zeros before Fourier transforming to obtain a smooth spectrum (reduce frequency step).
            Rt_padded = np.pad(Rt * damping, ((0, npad), (0, npad)))
            # Sine transform over first time axis (sine trasnform is equivalent to the imaginary part of the Fourier transform).
            St_partial = np.imag(fftpack.fft(Rt_padded, axis=0))
            # Sine transform over the second time axis.
            St_raw = np.imag(fftpack.fft(St_partial, axis=1))
            # Take the part of the spectrum that we want to plot, multiply by dt**2 to get the right units and scaling. We divide by 10**6 for convenience here, so the units will be [10**6 a.u.].
            St = St_raw[: nplot + 1, : nplot + 1] * dt**2 / 10**6
            St = St / np.max(abs(St))
            plt.rc("text", usetex=True)
            plt.rc("font", **{"family": "serif", "serif": ["Computer Modern"], "size": 8})
            plt.rcParams["contour.negative_linestyle"] = "solid"
            plt.rcParams["contour.linewidth"] = 0.01

            vmin = -1
            vmax =  1
            levels = np.arange(vmin, vmax, 0.05)

            clp.plotone([], [], axes[j,i], xlabel=r"$\omega_1 / 2 \pi c\ [\rm{cm}^{-1}]$" if j == 4 else None, ylabel=r"$\omega_2 / 2 \pi c\ [\rm{cm}^{-1}]$" if i == 0 else None, showlegend=False)
            axes[j,i].text(0.03, 0.3, text_list[j]+"\n"+labels[i]+dipole_list[j], transform=axes[j,i].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
            axes[j,i].text(0.97, 0.97, idx[j][i], transform=axes[j,i].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")

            cs = axes[j,i].contourf(
                St.T,
                levels=levels,
                vmin=vmin,
                vmax=vmax,
                cmap=mpl.cm.bwr,
                extend="both",
                extent=(wmin, wmax, wmin, wmax),
            )  # plots color map

    cax = fig.add_axes([0.91, 0.115, 1e-2, 0.76])  # [left, botton, width, heigth]
    cbar = fig.colorbar(cs, cax=cax)
    cbar.set_label('spectrum intensity [arb. units]', fontsize=12)
    cbar.set_ticks([-1, 0, 1])
    clp.adjust(savefile=f'./fig/2d-ir-raman.png')

def plot_2D_qc_Raman():
    
    fig, axes = clp.initialize(2, 2, width=4.3*2, height=4.3*0.85*2, LaTeX=True, fontsize=12, sharey=True, sharex=True, return_fig_args=True)
    # Useful constants.
    fstoau = 41.34137  # Femtosecond to a.u.
    autocm1 = 219474.63  # a.u. (Hartree energy unit) to wavenumber
    invfstocm1 = autocm1 / fstoau  # fs^-1 to cm^-1 conversion
    # Processing parameters
    dt = 1 * fstoau  # Timestep 1 fs
    nsteps = 250  # Number of steps in Rt
    npad = 1125*2  # Padding for FFT
    nplot = 175*2  # Number of frequency steps to plot.
    dw = 2 * np.pi / (nsteps + npad + 1) / dt * autocm1  # Frequency step for fft
    wmin = 0  # Minimum frequency in the plot
    wmax = nplot * dw  # Maximum frequency in the plot
    # Damping to avoid a hard cut off at long t_1/t_2.
    tau = 1e27
    damping = np.array(
        [
            [np.exp(-(((i + j) * dt / fstoau) ** 12) / tau) for i in range(nsteps + 1)]
            for j in range(nsteps + 1)
        ]
    )
    nlist = [64] 
    labels = [r"$N_{\rm{simu}}=%d$"%ni for ni in nlist]
    epsilon_list = ["0e-4","7e-4_did","7e-4_pc","7e-4_qc_l10ps"]
    text_list = ["outside-cavity",r"$\widetilde{\varepsilon}=7 \times 10^{-4}$ a.u.",r"$\widetilde{\varepsilon}=7 \times 10^{-4}$ a.u.",r"$\widetilde{\varepsilon}=7 \times 10^{-4}$ a.u."]
    idx = {0:["(a)","(b)"],1:["(c)","(d)"]}
    dipole_list = ["", r"DID-DID, $R_{\rm{m}}^{\rm{MD}}$", r"PC-DID, $R_{\rm{m}}^{\rm{MD}}$", r"DID-$\mathbf{q}_{\rm{c}}$, $R_{\rm{c}}^{\rm{MD}}$, $\tau=10$ ps"]
    
    for i in range(4):
        # Import Rt.
        Rt = np.loadtxt(f"./final_results/2d_iir_qc/2d-ir-raman_{epsilon_list[i]}.out")[:nsteps+1,:nsteps+1]
        # Rt = np.gradient(RtRaw,dt,axis=0,edge_order=2)
        # Pad with zeros before Fourier transforming to obtain a smooth spectrum (reduce frequency step).
        Rt_padded = np.pad(Rt * damping, ((0, npad), (0, npad)))
        # Sine transform over first time axis (sine trasnform is equivalent to the imaginary part of the Fourier transform).
        St_partial = np.imag(fftpack.fft(Rt_padded, axis=0))
        # Sine transform over the second time axis.
        St_raw = np.imag(fftpack.fft(St_partial, axis=1))
        # Take the part of the spectrum that we want to plot, multiply by dt**2 to get the right units and scaling. We divide by 10**6 for convenience here, so the units will be [10**6 a.u.].
        St = St_raw[: nplot + 1, : nplot + 1] * dt**2 / 10**6
        St = St / np.max(abs(St))
        plt.rc("text", usetex=True)
        plt.rc("font", **{"family": "serif", "serif": ["Computer Modern"], "size": 8})
        plt.rcParams["contour.negative_linestyle"] = "solid"
        plt.rcParams["contour.linewidth"] = 0.01

        vmin = -1
        vmax =  1
        levels = np.arange(vmin, vmax, 0.05)
        x0, y0 = i // 2, i % 2
        clp.plotone([], [], axes[x0, y0], xlabel=r"$\omega_1 / 2 \pi c\ [\rm{cm}^{-1}]$" if x0 == 1 else None, ylabel=r"$\omega_2 / 2 \pi c\ [\rm{cm}^{-1}]$" if y0 == 0 else None, showlegend=False)
        axes[x0, y0].text(0.03, 0.3, dipole_list[i]+"\n"+text_list[i]+"\n"+labels[0], transform=axes[x0, y0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
        axes[x0, y0].text(0.97, 0.97, idx[x0][y0], transform=axes[x0, y0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
    
        cs = axes[x0, y0].contourf(
            St.T,
            levels=levels,
            vmin=vmin,
            vmax=vmax,
            cmap=mpl.cm.bwr,
            extend="both",
            extent=(wmin, wmax, wmin, wmax),
        )  # plots color map

    for j in range(4):
        x0, y0 = j // 2, j % 2
        rect = mpatches.Rectangle((0,3000),1300,1000, fill=False,color="c",linewidth=2, linestyle='--')
        axes[x0, y0].add_patch(rect)        
  
    rect = mpatches.Rectangle((2740,2500),840*2,1650, fill=False,color="g",linewidth=2, linestyle='--')
    axes[1,0].add_patch(rect)
    rect = mpatches.Rectangle((3400,3000),600,900, fill=False,color="m",linewidth=2, linestyle='--')
    axes[1,0].add_patch(rect) 
    rect = mpatches.Rectangle((2800,2500),840*2,1650, fill=False,color="g",linewidth=2, linestyle='--')
    axes[0,1].add_patch(rect)
    rect = mpatches.Rectangle((3400,3000),400,900, fill=False,color="m",linewidth=2, linestyle='--')
    axes[0,1].add_patch(rect)   
    rect = mpatches.Rectangle((2740,2500),840*2,1650, fill=False,color="g",linewidth=2, linestyle='--')
    axes[1,1].add_patch(rect)

    cax = fig.add_axes([0.92, 0.115, 2e-2, 0.76])  # [left, botton, width, heigth]
    cbar = fig.colorbar(cs, cax=cax)
    cbar.set_label('spectrum intensity [arb. units]', fontsize=12)
    cbar.set_ticks([-1, 0, 1])
    clp.adjust(savefile=f'./fig/2d_qc_raman.png')

if __name__ == "__main__" :
    plot_qc()
    plot_IR()
    plot_IR_size()
    plot_Raman()
    plot_2D_IR_Raman()
    plot_2D_qc_Raman()