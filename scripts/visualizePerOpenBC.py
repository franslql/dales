import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pycwt as wavelet
from matplotlib import ticker, cm, colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

# Input
sigmaX = 0
sigmaT = 0
iz = 5
limitxy = 0.5
limitxz = 0.5
levels_thlmeanxz = np.linspace(-0.04,0.1,22)
levels_tkexz = np.linspace(0,1,25)
conf = 5
perc = 99
trun = 6*3600
category = 'u3r2'
pathPeriodic = f"../cases/periodic/"
pathOpenBC = f"../cases/openBC/x{sigmaX:03}y{sigmaX:03}z000t{sigmaT:03}/"
r = 2
# Set all fontsizes
SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18
plt.rc('font', size=MEDIUM_SIZE)         # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# -------------------- Plot crossxy ----------------------- #
# Load periodic data
nc   = Dataset(f"{pathPeriodic}fielddump/thl.000.nc",'r')
time = nc.variables['time'][:]; dt = time[1]-time[0]
time = time-time[0]+dt; time=time[time<=trun]; Nt = len(time); it=Nt-1
xt   = nc.variables['xt'][:]; Nx = len(xt); dx = xt[1]-xt[0]
yt   = nc.variables['yt'][:]; Ny = len(yt); dy = yt[1]-yt[0]
zt   = nc.variables['zt'][:]; Nz = len(zt); dz = zt[1]-zt[0]
thlxy_per = nc.variables['thl'][it,iz,:,:]; thlxy_mean = np.mean(thlxy_per)
thlprof_per = np.mean(nc.variables['thl'][it,:,:,:],axis=(1,2))
nc.close()
nc   = Dataset(f"{pathPeriodic}fielddump/u.000.nc",'r')
uprof_per = np.mean(nc.variables['u'][it,:,:,:],axis=(1,2))
nc.close()
nc   = Dataset(f"{pathPeriodic}tmser.000.nc",'r')
zi   = nc.variables['zi'][it]
nc.close()
# Load openBC data
nc   = Dataset(f"{pathOpenBC}fielddump/thl.001.nc",'r')
thlxy = nc.variables['thl'][it,iz,:,:]
thlprof = np.mean(nc.variables['thl'][it,:,:,:],axis=(1,2))
nc.close()
nc   = Dataset(f"{pathOpenBC}fielddump/u.001.nc",'r')
uprof = np.mean(nc.variables['u'][it,:,:,:],axis=(1,2))
nc.close()
# Plot crossxy
fig = plt.figure()
gs = fig.add_gridspec(2,3)
axs = [None]*3
axs[0] = fig.add_subplot(gs[0, 1:])
axs[1] = fig.add_subplot(gs[1, 1:])
axs[2] = fig.add_subplot(gs[:, 0])
#axs = gs.subplots()
C= axs[0].contourf(xt/1000,yt/1000,thlxy_per-thlxy_mean,levels=np.linspace(-limitxy,limitxy,25),extend='both',cmap='coolwarm')
axs[1].contourf(xt/1000,yt/1000,thlxy-thlxy_mean,levels=np.linspace(-limitxy,limitxy,25),extend='both',cmap='coolwarm')
#fig.suptitle(f"time = {np.round(time[it]/3600,1)} hours, z = {zt[iz]} m",fontsize=18)
axs[0].set_title('Periodic')
axs[1].set_title('OpenBC')
axs[0].set_xticks([])
axs[1].set_xticks([xt.min()/1000,3,6,9,12,15])
axs[1].set_xticklabels([0,3,6,9,12,15])
axs[0].set_yticks([1,2,3])
axs[1].set_yticks([1,2,3])
axs[1].set_xlabel('x (km)')
axs[1].text(-0.05*xt[-1]/1000,yt[-1]*1.13/1000,'y (km)',rotation='vertical',va='center',ha='center')
cb_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
cbar = fig.colorbar(C, cax=cb_ax,ticks=np.linspace(-limitxy,limitxy,5))
cbar.set_label(f"$\\theta-\\left<\\theta_{{per}}\\right>$ $(K)$")
#axs[2].plot(thlprof_per,zt,'k')
axs[2].plot(thlprof-thlprof_per,zt,label=r'$\left<\theta\right>-\left<\theta_{per}\right>$ (K)')
axs[2].plot(uprof-uprof_per,zt,label=r'$\left<u\right>-\left<u_{per}\right>$ (m/s)')
axs[2].plot([-.2,.2],[zi,zi],'--',color='0.6',label='Inversion layer')
axs[2].set_ylabel('z (m)')
axs[2].set_yticks([300,600,900,1200,1500,1800])
axs[2].set_xticks([-0.2,-0.1,0,0.1,0.2])
axs[2].legend(loc='upper left')
fig.set_size_inches(16, 9)
fig.savefig(f"../results/periodic_vs_openBC/thlxy.png",
           bbox_inches='tight', dpi=200)
# -------------------- Plot crossxz ----------------------- #
# Load periodic data
nc   = Dataset(f"{pathPeriodic}fielddump/thl.000.nc",'r')
iy   = int(Ny/2)
thlxz_per = nc.variables['thl'][it,:,iy,:]
thlxz_mean = np.mean(nc.variables['thl'][it,:,:,:],axis=(1,2))
nc.close()
# Load openBC data
nc   = Dataset(f"{pathOpenBC}fielddump/thl.001.nc",'r')
thlxz = nc.variables['thl'][it,:,iy,:]
nc.close()
# Plot crossxy
fig = plt.figure()
gs = fig.add_gridspec(2,1)
axs = gs.subplots()
C= axs[0].contourf(xt/1000,zt,thlxz_per-thlxz_mean[:,None],levels=np.linspace(-limitxz,limitxz,25),extend='both',cmap='coolwarm')
axs[1].contourf(xt/1000,zt,thlxz-thlxz_mean[:,None],levels=np.linspace(-limitxz,limitxz,25),extend='both',cmap='coolwarm')
axs[1].plot(xt[xt/r<1000]/1000,xt[xt/r<1000]/r,'k',alpha=0.5)
#fig.suptitle(f"time = {np.round(time[it]/3600,1)} hours, y = {np.round(yt[iy]/1000)} km",fontsize=18)
axs[0].set_title('Periodic')
axs[1].set_title('OpenBC')
axs[0].set_xticks([])
axs[1].set_xticks([xt.min()/1000,3,6,9,12,15])
axs[1].set_xticklabels([0,3,6,9,12,15])
axs[0].set_yticks([300,600,900,1200,1500,1800])
axs[1].set_yticks([300,600,900,1200,1500,1800])
axs[1].set_xlabel('x (km)')
axs[1].text(-0.06*xt[-1]/1000,zt[-1]*1.13,'z (m)',rotation='vertical',va='center',ha='center')
cb_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
cbar = fig.colorbar(C, cax=cb_ax)
cbar.set_label(f"$\\theta-\\left<\\theta_{{per}}\\right>$ $(K)$")
fig.set_size_inches(16, 9)
fig.savefig(f"../results/periodic_vs_openBC/thlxz.png",
           bbox_inches='tight', dpi=200)
# -------------------- Plot meancrossxz ----------------------- #
# Load periodic data
nc   = Dataset(f"{pathPeriodic}meancrossxz/thlwxzmean.000.nc",'r')
thlwxz_per = nc.variables['thlwxz'][it,:,:]
thlwxz_perc_high = np.percentile(thlwxz_per,100-conf/2,axis=1)
thlwxz_perc_low  = np.percentile(thlwxz_per,conf/2,axis=1)
zm = nc.variables['zm'][:]
nc.close()
nc   = Dataset(f"{pathPeriodic}meancrossxz/u2xzmean.000.nc",'r')
u2xz_per = nc.variables['u2xz'][it,:,:]
u2xz_per = np.append(u2xz_per,u2xz_per[:,0,None],axis=1)
nc.close()
nc   = Dataset(f"{pathPeriodic}meancrossxz/v2xzmean.000.nc",'r')
v2xz_per = nc.variables['v2xz'][it,:,:]
nc.close()
nc   = Dataset(f"{pathPeriodic}meancrossxz/w2xzmean.000.nc",'r')
w2xz_per = nc.variables['w2xz'][it,:,:]
w2xz_per[-1,:]=0.
nc.close()
# Load openBC data
nc   = Dataset(f"{pathOpenBC}meancrossxz/thlwxzmean.001.nc",'r')
thlwxz = nc.variables['thlwxz'][it,:,:]
nc.close()
nc   = Dataset(f"{pathOpenBC}meancrossxz/u2xzmean.001.nc",'r')
u2xz = nc.variables['u2xz'][it,:,:]
xm = nc.variables['xm'][:]
nc.close()
nc   = Dataset(f"{pathOpenBC}meancrossxz/v2xzmean.001.nc",'r')
v2xz = nc.variables['v2xz'][it,:,:]
nc.close()
nc   = Dataset(f"{pathOpenBC}meancrossxz/w2xzmean.001.nc",'r')
w2xz = nc.variables['w2xz'][it,:,:]
nc.close()
tke_per = 0.5*(0.5*(u2xz_per[:,:-1]+u2xz_per[:,1:])+v2xz_per+0.5*(w2xz_per[:-1,:]+w2xz_per[1:,:]))
tke_perc_high = np.percentile(tke_per,100-conf/2,axis=1)
tke_perc_low  = np.percentile(tke_per,conf/2,axis=1)
tke     = 0.5*(0.5*(u2xz[:,:-1]+u2xz[:,1:])+v2xz+0.5*(w2xz[:-1,:]+w2xz[1:,:]))
# Plot meancrossxz
fig = plt.figure()
gs = fig.add_gridspec(2,1)
axs = gs.subplots()
C= axs[0].contourf(xt/1000,zm,thlwxz_per,levels=levels_thlmeanxz,extend='both',cmap='viridis')
axs[0].contour(xt/1000,zm,thlwxz_per-thlwxz_perc_low[:,None],levels=[0],linestyles='dotted',colors='0.6')
axs[0].contour(xt/1000,zm,thlwxz_per-thlwxz_perc_high[:,None],levels=[0],linestyles='dashed',colors='0.6')
axs[1].contourf(xt/1000,zm,thlwxz,levels=levels_thlmeanxz,extend='both',cmap='viridis')
axs[1].contour(xt/1000,zm,thlwxz-thlwxz_perc_low[:,None],levels=[0],linestyles='dotted',colors='0.6')
axs[1].contour(xt/1000,zm,thlwxz-thlwxz_perc_high[:,None],levels=[0],linestyles='dashed',colors='0.6')
#axs[1].plot(xt[xt/r<1000]/1000,xt[xt/r<1000]/r,'k',alpha=0.5)
#fig.suptitle(f"time = {np.round(time[it]/3600,1)} hours",fontsize=18)
axs[0].set_title('Periodic')
axs[1].set_title('OpenBC')
axs[0].set_xticks([])
axs[1].set_xticks([xt.min()/1000,3,6,9,12,15])
axs[1].set_xticklabels([0,3,6,9,12,15])
axs[0].set_yticks([300,600,900,1200,1500,1800])
axs[1].set_yticks([300,600,900,1200,1500,1800])
axs[1].set_xlabel('x (km)')
axs[1].text(-0.06*xt[-1]/1000,zt[-1]*1.13,'z (m)',rotation='vertical',va='center',ha='center')
cb_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
cbar = fig.colorbar(C, cax=cb_ax)
cbar.set_label(f"$\\overline{{w'\\theta'}}$  (K m/s)")
fig.set_size_inches(16, 9)
fig.savefig(f"../results/periodic_vs_openBC/thlwxz.png",
           bbox_inches='tight', dpi=200)

fig = plt.figure()
gs = fig.add_gridspec(2,1)
axs = gs.subplots()
C= axs[0].contourf(xt/1000,zt,tke_per,levels=levels_tkexz,extend='both',cmap='viridis')
axs[0].contour(xt/1000,zt,tke_per-tke_perc_low[:,None],levels=[0],linestyles='dotted',colors='0.6')
axs[0].contour(xt/1000,zt,tke_per-tke_perc_high[:,None],levels=[0],linestyles='dashed',colors='0.6')
axs[1].contourf(xt/1000,zt,tke,levels=levels_tkexz,extend='both',cmap='viridis')
axs[1].contour(xt/1000,zt,tke-tke_perc_low[:,None],levels=[0],linestyles='dotted',colors='0.6')
axs[1].contour(xt/1000,zt,tke-tke_perc_high[:,None],levels=[0],linestyles='dashed',colors='0.6')
axs[1].plot(xt[xt/r<1000]/1000,xt[xt/r<1000]/r,'k',alpha=0.5)
#fig.suptitle(f"time = {np.round(time[it]/3600,1)} hours",fontsize=18)
axs[0].set_title('Periodic')
axs[1].set_title('OpenBC')
axs[0].set_xticks([])
axs[1].set_xticks([xt.min()/1000,3,6,9,12,15])
axs[1].set_xticklabels([0,3,6,9,12,15])
axs[0].set_yticks([300,600,900,1200,1500,1800])
axs[1].set_yticks([300,600,900,1200,1500,1800])
axs[1].set_xlabel('x (km)')
axs[1].text(-0.06*xt[-1]/1000,zt[-1]*1.13,'z (m)',rotation='vertical',va='center',ha='center')
cb_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
cbar = fig.colorbar(C, cax=cb_ax,ticks=np.linspace(0,1,5))
cbar.set_label(f"tke ($m^2/s^2$)")
fig.set_size_inches(16, 9)
fig.savefig(f"../results/periodic_vs_openBC/tkexz.png",
           bbox_inches='tight', dpi=200)
# -------------------- Plot wavelet results ----------------------- #
# Wavelet analysis
thlxy_diff_per = thlxy_per-np.mean(thlxy_per,axis=1)[:,None]
thlxy_diff     = thlxy-np.mean(thlxy,axis=1)[:,None]
# Normalize
thlxy_diff_per = thlxy_diff_per/np.mean(np.std(thlxy_diff_per,axis=1))
thlxy_diff = thlxy_diff/np.mean(np.std(thlxy_diff,axis=1))
# Set up continuous wavelet transform input
mother = wavelet.Morlet(4)
s0 = 2 * dx  # 2*resolution (nyquist frequency)
dj = 1 / 32  # Thirty two sub-octaves per octaves
J = int(5 / dj)  # Five powers of two with dj sub-octaves
# Calculate wavelet transform for every line
power_per = np.zeros((J+1,Nx))
power     = np.zeros((J+1,Nx))
alpha = 0
# fft_power = np.zeros(int(Nx/2)-1)
for iy in range(Ny):
    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(thlxy_diff_per[iy,:], dx, dj, s0, J, mother)
    power_per = power_per + np.abs(wave)**2/Ny
    alpha_temp, _, _ = wavelet.ar1(thlxy_diff_per[iy,:])  # Lag-1 autocorrelation for red noise
    alpha = alpha + alpha_temp/Ny
for iy in range(Ny):
    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(thlxy_diff[iy,:], dx, dj, s0, J, mother)
    power = power + np.abs(wave)**2/Ny
    # fft_power = fft_power + np.abs(fft)**2/Ny
period = 1/freqs
# Adjust for bias
power_per /= scales[:, None]
power /= scales[:, None]
# # Power spectra significance test
# signif, fft_theor = wavelet.significance(1.0, dx, scales, 0, alpha,
#                                          significance_level=0.01,
#                                          wavelet=mother)
# sig95_per = np.ones([1, Nx]) * signif[:, None]
# sig95_per = power_per / sig95_per
# sig95 = np.ones([1, Nx]) * signif[:, None]
# sig95 = power / sig95

power_perc_high = np.percentile(power_per,100-conf/2,axis=1)
power_perc_low  = np.percentile(power_per,conf/2,axis=1)
# Create levels colorbar
basepower = -4
npowers   = 3
levs    = np.zeros(npowers*9+1)
for ipower in range(npowers):
    levs[ipower*9:(ipower+1)*9] = np.arange(1,10)*10**(basepower+ipower)
levs[-1] = 10**(basepower+npowers)
# Create colormap
skip = 2
Greys = cm.get_cmap('Greys', 10+skip)
Purples = cm.get_cmap('Purples', 10+skip)
Blues = cm.get_cmap('Blues', 10+skip)
Greens = cm.get_cmap('Greens', 10+skip)
Oranges = cm.get_cmap('Oranges', 10+skip)
Reds = cm.get_cmap('Reds', 10+skip)
cmap = np.concatenate([#Greys(range(10+skip))[skip:,:],
                       #Purples(range(10+skip))[skip:,:],
                       Blues(range(10+skip))[skip:,:],
                       Greens(range(10+skip))[skip:,:],
                       #Oranges(range(10+skip))[skip:,:],
                       Reds(range(10+skip))[skip:,:]
                       ],axis=0)
cmap_in = ListedColormap(cmap)
cbarticklabels= [None]*(int(npowers*9)+1)
for ilabel in range(len(cbarticklabels)):
    if(ilabel%9==0):
        cbarticklabels[ilabel] = f"$10^{{{int(basepower+ilabel/9)}}}$"
    else:
        cbarticklabels[ilabel] = ''
# Plot results
fig = plt.figure()
gs = fig.add_gridspec(2,1)
axs = gs.subplots()
norm = colors.BoundaryNorm(np.log10(levs), cmap_in.N)
extent = [xt.min(),xt.max(),period.min(),period.max()]
C=axs[0].contourf(xt/1000,np.log10(period),np.log10(power_per),np.log10(levs),cmap=cmap_in,extend='both', norm=norm)
axs[0].contour(xt/1000,np.log10(period),np.log10(power_per)-np.log10(power_perc_low[:,None]),levels=[0],linestyles='dotted',colors='0.6',extent=extent)
axs[0].contour(xt/1000,np.log10(period),np.log10(power_per)-np.log10(power_perc_high[:,None]),levels=[0],linestyles='dashed',colors='0.6',extent=extent)
# axs[0].contour(xt/1000, np.log10(period), sig95_per, [-99, 1], colors='k', linewidths=2,
#            extent=[xt.min(), xt.max(), 0, max(period)])
axs[0].fill(np.concatenate([xt, xt[-1:] + dx, xt[-1:] + dx,
                           xt[:1] - dx, xt[:1] - dx])/1000,
        np.log10(np.concatenate([coi, [2*dx], period[-1:],
                           period[-1:], [2*dx]])),
        'k', alpha=0.3, hatch='x')
axs[0].set_xlim([xt.min()/1000,xt[-1]/1000])
axs[0].set_ylim(np.log10([period.min(),period.max()]))
axs[1].contourf(xt/1000,np.log10(period),np.log10(power),np.log10(levs),cmap=cmap_in,extend='both', norm=norm)
# axs[1].contour(xt/1000, np.log10(period), sig95, [-99, 1], colors='k', linewidths=2,
#            extent=[xt.min(), xt.max(), 0, max(period)])
axs[1].contour(xt/1000,np.log10(period),np.log10(power)-np.log10(power_perc_low[:,None]),levels=[0],linestyles='dotted',colors='0.6',extent=extent)
axs[1].contour(xt/1000,np.log10(period),np.log10(power)-np.log10(power_perc_high[:,None]),levels=[0],linestyles='dashed',colors='0.6',extent=extent)
axs[1].fill(np.concatenate([xt, xt[-1:] + dx, xt[-1:] + dx,
                           xt[:1] - dx, xt[:1] - dx])/1000,
        np.log10(np.concatenate([coi, [2*dx], period[-1:],
                           period[-1:], [2*dx]])),
        'k', alpha=0.3, hatch='x')
axs[1].set_xlim([xt.min()/1000,xt[-1]/1000])
axs[1].set_ylim(np.log10([period.min(),period.max()]))
# Set set_yticklabels
basepower = 2
npowers   = 2
skipStart = 1
skipEnd   = 6
yticks    = np.zeros(npowers*9+1)
yticklabels= [None]*(int(npowers*9)+1)
for ipower in range(npowers):
    yticks[ipower*9:(ipower+1)*9] = np.arange(1,10)*10**(basepower+ipower)
yticks[-1] = 10**(basepower+npowers)
for ilabel in range(len(yticklabels)):
    if(ilabel%9==0):
        yticklabels[ilabel] = f"$1\cdot 10^{int(basepower+ilabel/9)}$"
    else:
        yticklabels[ilabel] = ''
yticks = yticks[skipStart:-skipEnd]
yticklabels = yticklabels[skipStart:-skipEnd]
yticklabels[0] = f"${skipStart+1}\cdot 10^{basepower}$"
yticklabels[-1] = f"${10-skipEnd}\cdot 10^{basepower+npowers-1}$"
axs[0].set_yticks(np.log10(yticks))
axs[0].set_yticklabels(yticklabels)
axs[1].set_yticks(np.log10(yticks))
axs[1].set_yticklabels(yticklabels)
#fig.suptitle(f"time = {np.round(time[it]/3600,1)} hours, z = {zt[iz]} m",fontsize=18)
axs[0].set_title('Periodic')
axs[1].set_title('OpenBC')
axs[0].set_xticks([])
axs[1].set_xticks([xt.min()/1000,3,6,9,12,15])
axs[1].set_xticklabels([0,3,6,9,12,15])
axs[1].set_xlabel('x (km)')
axs[1].text(-0.08*xt[-1]/1000,3.9,r'$\lambda (m)$',rotation='vertical',va='center',ha='center')
cb_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
cbar = fig.colorbar(
    cm.ScalarMappable(cmap=cmap_in, norm=norm),
    cax=cb_ax,
    boundaries=np.log10(np.concatenate([[9*10**(-6),],levs,[2*10**(0),]])),
    extend='both',
    ticks=np.log10(levs),
    spacing='proportional')
cbar.ax.set_yticklabels(cbarticklabels)
cbar.set_label(r'$P_{\theta}$')
fig.set_size_inches(16, 9)
fig.savefig(f"../results/periodic_vs_openBC/thlWavelet.png",
           bbox_inches='tight', dpi=200)
