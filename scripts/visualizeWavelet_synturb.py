import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pycwt as wavelet
from matplotlib import ticker, cm, colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

# Input
sigmaX_array = np.array([0, 2, 4, 8, 16])
sigmaT_array = np.array([0, 6, 30, 180])
trun = 6*3600
conf = 5
iz = 5
N1 = len(sigmaX_array)
N2 = len(sigmaT_array)
pathCases = "../cases/"
pathPeriodic   = f"{pathCases}periodic/fielddump/"
pathOpenBC     = f"{pathCases}openBC/"
pathOpenBC_synturb     = f"{pathCases}openBC_synturb/"
# Set all fontsizes
SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18
plt.rc('font', size=MEDIUM_SIZE)         # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
# Load periodic data
nc   = Dataset(f"{pathPeriodic}thl.000.nc",'r')
time = nc.variables['time'][:]; dt = time[1]-time[0]
time = time-time[0]+dt; time=time[time<=trun]; Nt = len(time); it=Nt-1
xt   = nc.variables['xt'][:]; Nx = len(xt); dx = xt[1]-xt[0]; lx = (Nx+1)*dx
yt   = nc.variables['yt'][:]; Ny = len(yt); dy = yt[1]-yt[0]
zt   = nc.variables['zt'][:]; Nz = len(zt); dz = zt[1]-zt[0]
thl_per = nc.variables['thl'][it,iz,:,:]
nc.close()
# Set up continuous wavelet transform input
mother = wavelet.Morlet(4)
s0 = 2 * dx  # 2*resolution (nyquist frequency)
dj = 1 / 32  # Thirty two sub-octaves per octaves
J = int(5 / dj)  # Five powers of two with dj sub-octaves
# Periodic wavelet transform
# Normalize periodc data
thlDiff_per = thl_per-np.mean(thl_per,axis=1)[:,None]
thlDiff_per = thlDiff_per/np.mean(np.std(thlDiff_per,axis=1))
# Calculate wavelet transform for every line
power_per = np.zeros((J+1, Nx))
for iy in range(Ny):
    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(thlDiff_per[iy,:], dx, dj, s0, J, mother)
    power_per = power_per + np.abs(wave)**2/Ny
# Adjust for bias
power_per /= scales[:, None]
power_perc_high = np.percentile(power_per,100-conf/2,axis=1)
power_perc_low  = np.percentile(power_per,conf/2,axis=1)
# Do analysis for smoothed input simulations
fig = plt.figure()
gs = fig.add_gridspec(N1, N2, hspace=0.1, wspace=0.05)
axs = gs.subplots()
for i in range(N1):
    sigmaX = sigmaX_array[i]
    for j in range(N2):
        sigmaT = sigmaT_array[j]
        # Load and normalize data
        if(sigmaX == 0 and sigmaT == 0):
            pathData = f"{pathOpenBC}x{sigmaX:03}y{sigmaX:03}z000t{sigmaT:03}/fielddump/"
        else:
            pathData = f"{pathOpenBC_synturb}x{sigmaX:03}y{sigmaX:03}z000t{sigmaT:03}/fielddump/"
        nc = Dataset(
            f"{pathData}thl.001.nc", 'r')
        thl = nc.variables['thl'][it, iz, :, :]
        nc.close()
        # Normalize
        thlDiff = thl-np.mean(thl,axis=1)[:,None]
        thlDiff = thlDiff/np.mean(np.std(thlDiff,axis=1))
        # Calculate wavelet transform for every line
        power = np.zeros((J+1, Nx))
        for iy in range(np.shape(thlDiff)[0]):
            wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(thlDiff[iy, :], dx, dj, s0, J,
                                                                  mother)
            power = power + np.abs(wave)**2/Ny
        period = 1/freqs
        # Adjust for bias
        power /= scales[:, None]
        # Plot results
        # Create levels colorbar
        basepower = -4
        npowers = 3
        levs = np.zeros(npowers*9+1)
        for ipower in range(npowers):
            levs[ipower*9:(ipower+1)*9] = np.arange(1, 10) * \
                           10**(basepower+ipower)
        levs[-1] = 10**(basepower+npowers)
        # Create colormap
        skip = 2
        Greys = cm.get_cmap('Greys', 10+skip)
        Purples = cm.get_cmap('Purples', 10+skip)
        Blues = cm.get_cmap('Blues', 10+skip)
        Greens = cm.get_cmap('Greens', 10+skip)
        Oranges = cm.get_cmap('Oranges', 10+skip)
        Reds = cm.get_cmap('Reds', 10+skip)
        cmap = np.concatenate([#Greys(range(10+skip))[skip:, :],
                               #Purples(range(10+skip))[skip:, :],
                               Blues(range(10+skip))[skip:, :],
                               Greens(range(10+skip))[skip:, :],
                               #Oranges(range(10+skip))[skip:,:],
                               Reds(range(10+skip))[skip:, :]], axis=0)
        cmap_in = ListedColormap(cmap)
        cbarticklabels = [None]*(int(npowers*9)+1)
        for ilabel in range(len(cbarticklabels)):
            if(ilabel % 9 == 0):
                cbarticklabels[ilabel] = f"$10^{{{int(basepower+ilabel/9)}}}$"
            else:
                cbarticklabels[ilabel] = ''
        norm = colors.BoundaryNorm(np.log10(levs), cmap_in.N)
        C = axs[i, j].contourf(xt/1000, np.log10(period), np.log10(power),
                               np.log10(levs), cmap=cmap_in, extend='both', norm=norm)
        axs[i, j].contour(xt/1000, np.log10(period), np.log10(power)-np.log10(power_perc_low[:,None]),
                               levels=[0],linestyles='dotted',colors='0.4')
        axs[i, j].contour(xt/1000, np.log10(period), np.log10(power)-np.log10(power_perc_high[:,None]),
                               levels=[0],linestyles='dashed',colors='0.4')
        axs[i, j].fill(np.concatenate([xt, xt[-1:] + dx, xt[-1:] + dx,
                                      xt[:1] - dx, xt[:1] - dx])/1000,
                       np.log10(np.concatenate([coi, [2*dx], period[-1:],
                                               period[-1:], [2*dx]])),
                       'k', alpha=0.3, hatch='x')
        axs[i, j].set_xlim([0, xt[-1]/1000])
        axs[i, j].set_ylim(np.log10([period.min(), period.max()]))
        if(i < N1-1):
            axs[i, j].set_xticks([])
        else:
            axs[i,j].set_xticks([xt.min()/1000,3,6,9,12,15])
            axs[i,j].set_xticklabels([0,3,6,9,12,15])
        if(j < N2-1):
            axs[i, j].set_yticks([])
        else:
            basepower = 2
            npowers = 2
            skipStart = 1
            skipEnd = 6
            yticks = np.zeros(npowers*9+1)
            yticklabels = [None]*(int(npowers*9)+1)
            for ipower in range(npowers):
                yticks[ipower*9:(ipower+1)*9] = np.arange(1,
                                                          10)*10**(basepower+ipower)
            yticks[-1] = 10**(basepower+npowers)
            for ilabel in range(len(yticklabels)):
                if(ilabel % 9 == 0):
                    yticklabels[ilabel] = f"$1\cdot 10^{int(basepower+ilabel/9)}$"
                else:
                    yticklabels[ilabel] = ''
            yticks = yticks[skipStart:-skipEnd]
            yticklabels = yticklabels[skipStart:-skipEnd]
            if(skipStart != 0):
                yticklabels[0] = f"${skipStart+1}\cdot 10^{basepower}$"
            if(skipEnd != 0):
                yticklabels[-1] = f"${10-skipEnd}\cdot 10^{basepower+npowers-1}$"
            axs[i, j].set_yticks(np.log10(yticks))
            axs[i, j].set_yticklabels(yticklabels)
            axs[i, j].yaxis.tick_right()
        if(j == 0):
            axs[i, j].set_ylabel(str(sigmaX)+r'$\Delta x$')
        if(i == 0):
            axs[i, j].set_title(str(sigmaT)+r'$\Delta t$')
axs[0, 0].text((-lx*0.15)/1000, -0.3, r'$\sigma_x$',
               rotation='vertical', va='center', ha='center')
axs[0, 0].text(((N2/2)*lx+(N2-1)/2*0.05*lx)/1000, 4.2,
               r'$\sigma_t$', va='center', ha='center')
axs[N1-1, 0].text(((N2/2)*lx+(N2-1)/2*0.05*lx)/1000, 1.7,
                  'x (km)', va='center', ha='center')
axs[0, N2-1].text(xt[-1]/1000*1.3, -0.3, r'$\lambda (m)$', rotation='vertical', va='center', ha='center')
cb_ax = fig.add_axes([0.065, 0.15, 0.015, 0.7])
cbar = fig.colorbar(
    cm.ScalarMappable(cmap=cmap_in, norm=norm),
    cax=cb_ax,
    boundaries=np.log10(np.concatenate([[9*10**(-6), ], levs, [2*10**(0), ]])),
    extend='both',
    ticks=np.log10(levs),
    spacing='proportional')
cbar.ax.set_yticklabels(cbarticklabels)
cbar.set_label(r'$P_{\theta}$', size=14)
cb_ax.yaxis.set_ticks_position('left')
cb_ax.yaxis.set_label_position('left')
fig.set_size_inches(16, 9)
#fig.suptitle(
#    f"time = {int(np.round(time[it]/3600,0))} hours, z = {int(np.round(zt[iz],0))} m", fontsize=18)
plt.savefig(
    f"../results/openBC_synturb/thlWavelet.png",
               bbox_inches='tight', dpi=200)
