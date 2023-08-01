import numpy as np
import scipy.ndimage as sp
import funcDALES as DALES
from netCDF4 import Dataset

# Input
pathCases = "../cases/"
pathWrite = "../input/boundary_input/"
sigmaX_array = np.array([2])
sigmaT_array = np.array([0])
runtime = 21600
pathPeriodic = f"{pathCases}periodic/"
# Read periodic data
# profile input
data = np.loadtxt(f"{pathPeriodic}prof.inp.000")
zt   = data[:,0]
thl0 = data[:,1]
qt0  = data[:,2]
u0   = data[:,3]
v0   = data[:,4]
e120 = data[:,5]
ktot = len(zt)
# Read crossyz data
nc = Dataset(f"{pathPeriodic}crossyz/0002/thlyz.0002.000.nc",'r')
thlyz = nc.variables['thlyz'][:,:,:]
yt = nc.variables['yt'][:]
jtot = len(yt)
time = nc.variables['time'][:]
time = time-2*time[0]+time[1]
nt = len(time)
nc.close()
nc = Dataset(f"{pathPeriodic}crossyz/0002/qtyz.0002.000.nc",'r')
qtyz = nc.variables['qtyz'][:,:,:]
nc.close()
nc = Dataset(f"{pathPeriodic}crossyz/0002/uyz.0002.000.nc",'r')
uyz = nc.variables['uyz'][:,:,:]
nc.close()
nc = Dataset(f"{pathPeriodic}crossyz/0002/vyz.0002.000.nc",'r')
vyz = nc.variables['vyz'][:,:,:]
ym = nc.variables['ym'][:]
ym = np.append(ym,ym[-1]+ym[-1]-ym[-2])
nc.close()
nc = Dataset(f"{pathPeriodic}crossyz/0002/wyz.0002.000.nc",'r')
wyz = nc.variables['wyz'][:,:,:]
zm = nc.variables['zm'][:]
zm = np.append(zm,zm[-1]+(zm[-1]-zm[-2]))
nc.close()

# Read crossxz data
nc = Dataset(f"{pathPeriodic}crossxz/0002/thlxz.0002.000.nc",'r')
thlxz = nc.variables['thlxz'][:,:,:]
xt = nc.variables['xt'][:]
itot = len(xt)
nc.close()
nc = Dataset(f"{pathPeriodic}crossxz/0002/qtxz.0002.000.nc",'r')
qtxz = nc.variables['qtxz'][:,:,:]
nc.close()
nc = Dataset(f"{pathPeriodic}crossxz/0002/uxz.0002.000.nc",'r')
uxz = nc.variables['uxz'][:,:,:]
xm = nc.variables['xm'][:]
xm = np.append(xm,xm[-1]+(xm[-1]-xm[-2]))
nc.close()
nc = Dataset(f"{pathPeriodic}crossxz/0002/vxz.0002.000.nc",'r')
vxz = nc.variables['vxz'][:,:,:]
nc.close()
nc = Dataset(f"{pathPeriodic}crossxz/0002/wxz.0002.000.nc",'r')
wxz = nc.variables['wxz'][:,:,:]
nc.close()

# Read crossxy data
nc = Dataset(f"{pathPeriodic}crossxy/{ktot:04}/thlxy.{ktot:04}.000.nc",'r')
thlxy = nc.variables['thlxy'][:,:,:]
nc.close()
nc = Dataset(f"{pathPeriodic}crossxy/{ktot:04}/qtxy.{ktot:04}.000.nc",'r')
qtxy = nc.variables['qtxy'][:,:,:]
nc.close()
nc = Dataset(f"{pathPeriodic}crossxy/{ktot:04}/uxy.{ktot:04}.000.nc",'r')
uxy = nc.variables['uxy'][:,:,:]
nc.close()
nc = Dataset(f"{pathPeriodic}crossxy/{ktot:04}/vxy.{ktot:04}.000.nc",'r')
vxy = nc.variables['vxy'][:,:,:]
nc.close()
wxy_filt = np.zeros((nt,jtot,itot))

grid = DALES.Grid(xm[-1], ym[-1], itot, jtot, ktot, zm[1]-zm[0], zm[1]-zm[0], 0., False)
# Filter periodic output and create open
for sigmaX in sigmaX_array:
    for sigmaT in sigmaT_array:
        sigmaY = sigmaX
        sigmaZ = 0
        experiment = f"x{sigmaX:03}y{sigmaY:03}z{sigmaZ:03}t{sigmaT:03}"
        print(f"Start creating input for experiment {experiment}")
        # Downscale spatiotemporal data using a gaussian smoother
        thlxz_filt = sp.gaussian_filter(thlxz,[sigmaT,sigmaZ,sigmaX],order=0,mode=['nearest','nearest','wrap'])
        qtxz_filt = sp.gaussian_filter(qtxz,[sigmaT,sigmaZ,sigmaX],order=0,mode=['nearest','nearest','wrap'])
        uxz_filt = sp.gaussian_filter(uxz,[sigmaT,sigmaZ,sigmaX],order=0,mode=['nearest','nearest','wrap'])
        vxz_filt = sp.gaussian_filter(vxz,[sigmaT,sigmaZ,sigmaX],order=0,mode=['nearest','nearest','wrap'])
        wxz_filt = sp.gaussian_filter(wxz,[sigmaT,sigmaZ,sigmaX],order=0,mode=['nearest','nearest','wrap'])

        thlyz_filt = sp.gaussian_filter(thlyz,[sigmaT,sigmaZ,sigmaY],order=0,mode=['nearest','nearest','wrap'])
        qtyz_filt = sp.gaussian_filter(qtyz,[sigmaT,sigmaZ,sigmaY],order=0,mode=['nearest','nearest','wrap'])
        uyz_filt = sp.gaussian_filter(uyz,[sigmaT,sigmaZ,sigmaY],order=0,mode=['nearest','nearest','wrap'])
        vyz_filt = sp.gaussian_filter(vyz,[sigmaT,sigmaZ,sigmaY],order=0,mode=['nearest','nearest','wrap'])
        wyz_filt = sp.gaussian_filter(wyz,[sigmaT,sigmaZ,sigmaY],order=0,mode=['nearest','nearest','wrap'])

        thlxy_filt = sp.gaussian_filter(thlxy,[sigmaT,sigmaY,sigmaX],order=0,mode=['nearest','wrap','wrap'])
        qtxy_filt = sp.gaussian_filter(qtxy,[sigmaT,sigmaY,sigmaX],order=0,mode=['nearest','wrap','wrap'])
        uxy_filt = sp.gaussian_filter(uxy,[sigmaT,sigmaY,sigmaX],order=0,mode=['nearest','wrap','wrap'])
        vxy_filt = sp.gaussian_filter(vxy,[sigmaT,sigmaY,sigmaX],order=0,mode=['nearest','wrap','wrap'])

        # Add extra half levels
        vyz_filt = np.append(vyz_filt,vyz_filt[:,:,0,None],axis=2)
        uxz_filt = np.append(uxz_filt,uxz_filt[:,:,0,None],axis=2)
        uxy_filt = np.append(uxy_filt,uxy_filt[:,:,0,None],axis=2)
        vxy_filt = np.append(vxy_filt,vxy_filt[:,0,None,:],axis=1)
        wxz_filt = np.append(wxz_filt,np.zeros((nt,1,itot)),axis=1)
        wyz_filt = np.append(wyz_filt,np.zeros((nt,1,jtot)),axis=1)
        wyz_filt[:,0,:] = 0.
        wxz_filt[:,0,:] = 0.

        # Create openboundaries.inp.iexnpr.nc
        openBC=DALES.OpenBoundariesInp(grid,iexpnr=experiment, path=pathWrite)
        # Write initial fields to file
        openBC.writeTime(0.,0)
        openBC.writeWest(thl0[:,None]*np.ones((ktot,jtot)),
                         qt0[:,None]*np.ones((ktot,jtot)),
                         u0[:,None]*np.ones((ktot,jtot)),
                         v0[:,None]*np.ones((ktot,jtot+1)),
                         np.zeros((ktot+1,jtot)),
                         e120[:,None]*np.ones((ktot,jtot)),0)
        openBC.writeEast(thl0[:,None]*np.ones((ktot,jtot)),
                         qt0[:,None]*np.ones((ktot,jtot)),
                         u0[:,None]*np.ones((ktot,jtot)),
                         v0[:,None]*np.ones((ktot,jtot+1)),
                         np.zeros((ktot+1,jtot)),
                         e120[:,None]*np.ones((ktot,jtot)),0)
        openBC.writeSouth(thl0[:,None]*np.ones((ktot,itot)),
                         qt0[:,None]*np.ones((ktot,itot)),
                         u0[:,None]*np.ones((ktot,itot+1)),
                         v0[:,None]*np.ones((ktot,itot)),
                         np.zeros((ktot+1,itot)),
                         e120[:,None]*np.ones((ktot,itot)),0)
        openBC.writeNorth(thl0[:,None]*np.ones((ktot,itot)),
                         qt0[:,None]*np.ones((ktot,itot)),
                         u0[:,None]*np.ones((ktot,itot+1)),
                         v0[:,None]*np.ones((ktot,itot)),
                         np.zeros((ktot+1,itot)),
                         e120[:,None]*np.ones((ktot,itot)),0)
        openBC.writeTop(thl0[-1]*np.ones((jtot,itot)),
                         qt0[-1]*np.ones((jtot,itot)),
                         u0[-1]*np.ones((jtot,itot+1)),
                         v0[-1]*np.ones((jtot+1,itot)),
                         np.zeros((jtot,itot)),
                         e120[-1]*np.ones((jtot,itot)),0)

        # Write later fields to file
        for it in range(nt):
            if(time[it]>runtime): break
            openBC.writeTime(time[it],it+1)
            openBC.writeWest(thlyz_filt[it,:,:],qtyz_filt[it,:,:],uyz_filt[it,:,:],vyz_filt[it,:,:],wyz_filt[it,:,:],e120[:,None]*np.ones((ktot,jtot)),it+1)
            openBC.writeEast(thlyz_filt[it,:,:],qtyz_filt[it,:,:],uyz_filt[it,:,:],vyz_filt[it,:,:],wyz_filt[it,:,:],e120[:,None]*np.ones((ktot,jtot)),it+1)
            openBC.writeSouth(thlxz_filt[it,:,:],qtxz_filt[it,:,:],uxz_filt[it,:,:],vxz_filt[it,:,:],wxz_filt[it,:,:],e120[:,None]*np.ones((ktot,itot)),it+1)
            openBC.writeNorth(thlxz_filt[it,:,:],qtxz_filt[it,:,:],uxz_filt[it,:,:],vxz_filt[it,:,:],wxz_filt[it,:,:],e120[:,None]*np.ones((ktot,itot)),it+1)
            openBC.writeTop(thlxy_filt[it,:,:],qtxy_filt[it,:,:],uxy_filt[it,:,:],vxy_filt[it,:,:],wxy_filt[it,:,:],np.zeros((jtot,itot)),it+1)
        # Close file
        openBC.exit()
        del openBC
