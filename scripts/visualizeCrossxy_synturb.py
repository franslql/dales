import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
# Input
sigmaX_array = np.array([0, 2, 4, 8, 16])
sigmaT_array = np.array([0, 6, 30, 180])
varname_array = ['thl','u','w']
label_array = [f"\\theta",'u','w']
units_array = ['K', 'm/s', 'm/s']
limit = np.array([0.5,2,2])
sigmaZ = 0
iz = 5
conf = 5
trun = 6*3600
r = 2
pathCases = "../cases/"
pathPeriodic   = f"{pathCases}periodic/fielddump/"
pathOpenBC_synturb    = f"{pathCases}openBC_synturb/"
pathOpenBC = f"{pathCases}openBC/" 
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
for ivar in range(len(varname_array)):
    varname = varname_array[ivar]
    units = units_array[ivar]
    label = label_array[ivar]
    nc = Dataset(f"{pathPeriodic}{varname}.000.nc",'r')
    time = nc.variables['time'][:]
    time = time-2*time[0]+time[1]
    it = np.where(time==trun)[0][0]
    varmean_per = np.mean(nc.variables[varname][it,iz,:,:])
    nc.close()
    fig = plt.figure()
    gs = fig.add_gridspec(len(sigmaX_array),len(sigmaT_array), hspace=0.1, wspace=0.05)#96/256*0.15)
    axs = gs.subplots()
    for i in range(len(sigmaX_array)):
        sigmaX = sigmaX_array[i]
        for j in range(len(sigmaT_array)):
            sigmaT = sigmaT_array[j]
            if(sigmaX==0 and sigmaT==0):
                pathData = f"{pathOpenBC}x{sigmaX:03}y{sigmaX:03}z{sigmaZ:03}t{sigmaT:03}/fielddump/"
            else:
                pathData = f"{pathOpenBC_synturb}x{sigmaX:03}y{sigmaX:03}z{sigmaZ:03}t{sigmaT:03}/fielddump/"
            nc = Dataset(f"{pathData}{varname}.001.nc",'r')
            time = nc.variables['time'][:]
            time = time-2*time[0]+time[1]
            it = np.where(time==trun)[0][0]
            var = nc.variables[varname][it,iz,:,:]
            for dimension in nc.dimensions.values():
                if(dimension.name[0]=='x'):
                    x = nc.variables[dimension.name][:]
                elif(dimension.name[0]=='y'):
                    y = nc.variables[dimension.name][:]
                elif(dimension.name[0]=='z'):
                    z = nc.variables[dimension.name][:]
            nc.close()
            # Plot crossection
            C=axs[i,j].contourf(x/1000,y/1000,var-varmean_per,levels=np.linspace(-limit[ivar],limit[ivar],25),extend='both',cmap='coolwarm')
            # Set axis labels and ticks
            if(i<len(sigmaX_array)-1):
                axs[i,j].set_xticks([])
            else:
                axs[i,j].set_xticks([x.min()/1000,3,6,9,12,15])
                axs[i,j].set_xticklabels([0,3,6,9,12,15])
            if(j<len(sigmaT_array)-1):
                axs[i,j].set_yticks([])
            else:
                axs[i,j].set_yticks([1,2,3])
                axs[i,j].yaxis.tick_right()
            if(j==0):
                axs[i,j].set_ylabel(str(sigmaX)+r'$\Delta x$')
            if(i==0):
                axs[i,j].set_title(str(sigmaT)+r'$\Delta t$')
    axs[0,0].text(-x[-1]*0.15/1000,-(len(sigmaX_array)/2-1)*y[-1]/1000-((len(sigmaX_array)-1)/2)*0.1*y[-1]/1000,r'$\sigma_x$',rotation='vertical',va='center',ha='center')
    axs[0,0].text((len(sigmaT_array)/2)*x[-1]/1000+(len(sigmaT_array)-1)/2*len(z)/len(x)*0.1*x[-1]/1000,y[-1]/1000+0.4*y[-1]/1000,r'$\sigma_t$',va='center',ha='center')
    axs[len(sigmaX_array)-1,0].text((len(sigmaT_array)/2)*x[-1]/1000+(len(sigmaT_array)-1)/2*len(z)/len(x)*0.1*x[-1]/1000,-y[-1]/1000*0.45,'x (km)',va='center',ha='center')
    axs[0,len(sigmaT_array)-1].text(x[-1]/1000*1.2,-(len(sigmaX_array)/2-1)*y[-1]/1000-((len(sigmaX_array)-1)/2)*0.1*y[-1]/1000,'y (km)',rotation='vertical',va='center',ha='center')
    cb_ax = fig.add_axes([0.065, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(C, cax=cb_ax,ticks=np.linspace(-limit[ivar],limit[ivar],5))
    cbar.set_label(f"${label}-\\left<{label}_{{per}}\\right>$ ({units})",size=14)
    cb_ax.yaxis.set_ticks_position('left')
    cb_ax.yaxis.set_label_position('left')
    fig.set_size_inches(16, 9)
    #fig.suptitle(f"time = {np.round(time[it]/3600,1)} hours, z = {np.round(z[iz],1)} m",fontsize=18)
    plt.savefig(f"../results/openBC_synturb/{varname}xy.png",
                bbox_inches='tight', dpi=200)
