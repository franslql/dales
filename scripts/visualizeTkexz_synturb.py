import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
# Input
trun = 6*3600
conf = 5
sigmaX_array = np.array([0, 2, 4, 8, 16])
sigmaT_array = np.array([0, 6, 30, 180])
sigmaZ = 0
pathCases = "../cases/"
pathPeriodic   = f"{pathCases}periodic/meancrossxz/"
pathOpenBC     = f"{pathCases}openBC/"
pathOpenBC_synturb = f"{pathCases}openBC_synturb/"
r = 2
levels_tkexz = np.linspace(0,1,25)
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
nc   = Dataset(f"{pathPeriodic}u2xzmean.000.nc",'r')
time = nc.variables['time'][:]; dt = time[1]-time[0]
time = time-time[0]+dt; time=time[time<=trun]; Nt = len(time); it=Nt-1
u2xz_per = nc.variables['u2xz'][it,:,:]
u2xz_per = np.append(u2xz_per,u2xz_per[:,0,None],axis=1)
nc.close()
nc   = Dataset(f"{pathPeriodic}/v2xzmean.000.nc",'r')
v2xz_per = nc.variables['v2xz'][it,:,:]
nc.close()
nc   = Dataset(f"{pathPeriodic}w2xzmean.000.nc",'r')
w2xz_per = nc.variables['w2xz'][it,:,:]
w2xz_per[-1,:]=0.
nc.close()
tke_per = 0.5*(0.5*(u2xz_per[:,:-1]+u2xz_per[:,1:])+v2xz_per+0.5*(w2xz_per[:-1,:]+w2xz_per[1:,:]))
tke_perc_high = np.percentile(tke_per,100-conf/2,axis=1)
tke_perc_low  = np.percentile(tke_per,conf/2,axis=1)
fig = plt.figure()
gs = fig.add_gridspec(len(sigmaX_array), len(
    sigmaT_array), hspace=0.1, wspace=0.05)#96/256*0.1)
axs = gs.subplots()
for i in range(len(sigmaX_array)):
    sigmaX = sigmaX_array[i]
    for j in range(len(sigmaT_array)):
        sigmaT = sigmaT_array[j]
        if(sigmaX == 0 and sigmaT == 0):
            pathData = f"{pathOpenBC}x{sigmaX:03}y{sigmaX:03}z{sigmaZ:03}t{sigmaT:03}/meancrossxz/"
        else:
            pathData = f"{pathOpenBC_synturb}x{sigmaX:03}y{sigmaX:03}z{sigmaZ:03}t{sigmaT:03}/meancrossxz/"
        nc = Dataset(f"{pathData}u2xzmean.001.nc", 'r')
        time = nc.variables['time'][:]
        time = time-2*time[0]+time[1]
        it = np.where(time==trun)[0][0]
        u2 = nc.variables['u2xz'][it,:,:]
        u2 = 0.5*(u2[:,:-1]+u2[:,1:])
        nc.close()
        nc = Dataset(f"{pathData}v2xzmean.001.nc", 'r')
        v2 = nc.variables['v2xz'][it,:,:]
        x = nc.variables['xt'][:]
        z = nc.variables['zt'][:]
        x = x + (x[1]-x[0])/2
        nc.close()
        nc = Dataset(f"{pathData}w2xzmean.001.nc", 'r')
        w2 = 0.5*(nc.variables['w2xz'][it,1:,:]+nc.variables['w2xz'][it,:-1,:])
        nc.close()
        tke = 0.5*(u2+v2+w2)
        # Plot tke
        C = axs[i,j].contourf(x/1000,z,tke,levels=levels_tkexz,extend='both',cmap='viridis')
        axs[i,j].contour(x/1000,z,tke-tke_perc_low[:,None],levels=[0],linestyles='dotted',colors='0.6')
        axs[i,j].contour(x/1000,z,tke-tke_perc_high[:,None],levels=[0],linestyles='dashed',colors='0.6')
        axs[i,j].plot(x[x/r<1000]/1000,x[x/r<1000]/r,'k',alpha=0.5)
        # Set axis labels and ticks
        if(i<len(sigmaX_array)-1):
            axs[i,j].set_xticks([])
        else:
            axs[i,j].set_xticks([x.min()/1000,3,6,9,12,15])
            axs[i,j].set_xticklabels([0,3,6,9,12,15])
        if(j<len(sigmaT_array)-1):
            axs[i,j].set_yticks([])
        else:
            axs[i,j].yaxis.tick_right()
            axs[i,j].set_yticks([300,900,1500])
        if(j==0):
            axs[i,j].set_ylabel(str(sigmaX)+r'$\Delta x$')
        if(i==0):
            axs[i,j].set_title(str(sigmaT)+r'$\Delta t$')
axs[0,0].text(-x[-1]*0.15/1000,-(len(sigmaX_array)/2-1)*z[-1]-((len(sigmaX_array)-1)/2)*0.1*z[-1],r'$\sigma_x$',rotation='vertical',va='center',ha='center')
axs[0,0].text((len(sigmaT_array)/2)*x[-1]/1000+(len(sigmaT_array)-1)/2*len(z)/len(x)*0.1*x[-1]/1000,z[-1]+0.4*z[-1],r'$\sigma_t$',va='center',ha='center')
axs[len(sigmaX_array)-1,0].text((len(sigmaT_array)/2)*x[-1]/1000+(len(sigmaT_array)-1)/2*len(z)/len(x)*0.1*x[-1]/1000,-z[-1]*0.45,'x (km)',va='center',ha='center')
axs[0,len(sigmaT_array)-1].text(x[-1]/1000*1.3,-(len(sigmaX_array)/2-1)*z[-1]-((len(sigmaX_array)-1)/2)*0.1*z[-1],'z (m)',rotation='vertical',va='center',ha='center')
cb_ax = fig.add_axes([0.065, 0.15, 0.015, 0.7])
cbar = fig.colorbar(C, cax=cb_ax,ticks=np.linspace(0,1,5))
cbar.set_label(f"tke ($m^2 s^{{-2}}$)",size=14)
cb_ax.yaxis.set_ticks_position('left')
cb_ax.yaxis.set_label_position('left')
fig.set_size_inches(16, 9)
#fig.suptitle(f"time = {np.round(time[it]/3600,1)} hours")
plt.savefig(f"../results/openBC_synturb/tkexz.png",
            bbox_inches='tight', dpi=200)
