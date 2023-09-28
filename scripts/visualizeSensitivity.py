import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt 
pathSensitivity = '../cases/sensitivity/'
pathPeriodic = '../cases/periodic/'
pathReference = '../cases/openBC/x000y000z000t000/'
tplot = 6*3600
tauh_array = np.array([0,60])
pbc_array  = np.array([2,4])
dint_array = np.array([0,5])
linestyles = ['--',':']
linewidth_ref = 2
linewidth = 1.5
markersize = 2.5
vars=['thl','u','wthlr','u2r']
titles=[f"$\\overline{{\\theta}}-\\overline{{\\theta}}_{{per}}$",\
        f"$\\overline{{u}}-\\overline{{u}}_{{per}}$",\
        f"$\\overline{{w'\\theta'}}-\\overline{{w'\\theta'}}_{{per}}$",\
        f"$\\overline{{u'^2}}-\\overline{{u'^2}}_{{per}}$"]
units=[f"$(K)$", f"$(ms^{{-1}})$", f"$(Kms^{{-1}})$", f"$(m^2s^{{-2}})$"]
limits = [[-1,1],[-0.5,0.5],[-0.03,0.03],[-0.7,0.7]]
legend_dint = [f"$\\Delta x_{{int}},\\Delta y_{{int}} = \\Delta x, \\Delta y$", f"$\\Delta x_{{int}},\\Delta y_{{int}} = 0.5L_x, 0.5L_y$"]
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
fig = plt.figure()
gs = fig.add_gridspec(25,4)
axs = [None]*4
axs[0] = fig.add_subplot(gs[:24,0])
axs[1] = fig.add_subplot(gs[:24,1])
axs[2] = fig.add_subplot(gs[:24,2])
axs[3] = fig.add_subplot(gs[:24,3])
for ivar in range(len(vars)):
  var = vars[ivar]
  # Open periodic simulation
  nc = Dataset(f"{pathPeriodic}profiles.000.nc",'r')
  z = nc.variables[nc.variables[var].dimensions[1]][:]
  time = nc.variables['time'][:]
  it = np.argwhere(time==tplot)[0][0]
  profPeriodic = nc.variables[var][it,:]
  nc.close()
  # Open integration length scale variation simulations
  for i in range(len(dint_array)):
    dint = dint_array[i] 
    nc = Dataset(f"{pathSensitivity}dint{dint:02}/profiles.001.nc",'r')
    z = nc.variables[nc.variables[var].dimensions[1]][:]
    time = nc.variables['time'][:]
    it = np.argwhere(time==tplot)[0][0]
    prof = nc.variables[var][it,:]
    axs[ivar].plot(prof-profPeriodic,z,f"{linestyles[i]}C0",markersize=markersize,linewidth=linewidth,label=legend_dint[i])
    nc.close()
  # Open pbc variation simulations
  for i in range(len(pbc_array)):
    pbc = pbc_array[i] 
    nc = Dataset(f"{pathSensitivity}pbc{pbc}/profiles.001.nc",'r')
    z = nc.variables[nc.variables[var].dimensions[1]][:]
    time = nc.variables['time'][:]
    it = np.argwhere(time==tplot)[0][0]
    prof = nc.variables[var][it,:]
    axs[ivar].plot(prof-profPeriodic,z,f"{linestyles[i]}C1",markersize=markersize,linewidth=linewidth,label=f"$p={pbc}$")
    nc.close()
  # Open tauh variation simulations
  for i in range(len(tauh_array)):
    tauh = tauh_array[i] 
    nc = Dataset(f"{pathSensitivity}tauh{tauh:03}/profiles.001.nc",'r')
    z = nc.variables[nc.variables[var].dimensions[1]][:]
    time = nc.variables['time'][:]
    it = np.argwhere(time==tplot)[0][0]
    prof = nc.variables[var][it,:]
    axs[ivar].plot(prof-profPeriodic,z,f"{linestyles[i]}C2",markersize=markersize,linewidth=linewidth,label=f"$\\tau_0 = {tauh}s$")
    nc.close()
  # Open simulation without buoyancy term at top boundary
  nc = Dataset(f"{pathSensitivity}lbuoyoff/profiles.001.nc",'r')
  z = nc.variables[nc.variables[var].dimensions[1]][:]
  time = nc.variables['time'][:]
  it = np.argwhere(time==tplot)[0][0]
  prof = nc.variables[var][it,:]
  axs[ivar].plot(prof-profPeriodic,z,f"{linestyles[0]}C3",linewidth=linewidth,label=f"Buoyancy top off")
  nc.close()
  # Open reference simulation
  nc = Dataset(f"{pathReference}profiles.001.nc",'r')
  zt = nc.variables[nc.variables[var].dimensions[1]][:]
  time = nc.variables['time'][:]
  it = np.argwhere(time==tplot)[0][0]
  prof = nc.variables[var][it,:]
  axs[ivar].plot(prof-profPeriodic,z,'k',linewidth=linewidth_ref,label=f"Default settings")
  nc.close()
  axs[ivar].grid()
  axs[ivar].set_xlim([limits[ivar][0],limits[ivar][1]])
  axs[ivar].set_ylim([0,2000])
  axs[ivar].set_xticks(np.linspace(limits[ivar][0],limits[ivar][1],5))
  for label in axs[ivar].xaxis.get_ticklabels()[1::2]:
    label.set_visible(False)
  if(ivar>0): 
    axs[ivar].set_yticklabels('')
  else:
    axs[ivar].set_ylabel('z (m)')
    h,l = axs[ivar].get_legend_handles_labels()
  axs[ivar].set_title(titles[ivar]+' '+units[ivar])
fig.legend(h,l,ncol=4,loc=8)
fig.set_size_inches(16, 9)
fig.savefig(f"../results/sensitivity.png",
           bbox_inches='tight', dpi=200)