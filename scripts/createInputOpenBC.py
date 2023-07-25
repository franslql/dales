import numpy as np
import scipy.ndimage as sp
import funcDALES as DALES
from netCDF4 import Dataset

# Input
# Experiment
pathStart = "../simulations/"
category = "u075r05"
sigmaX_array = np.array([4,8,16,32])
sigmaT_array = np.array([0,6,30,90,180,360])
runtime = 21600
pathWrite = f"{pathStart}{category}/openBC/input/openboundaries/"
pathPeriodic = f"{pathStart}/{category}/periodic/"
for sigmaX in sigmaX_array:
    for sigmaT in sigmaT_array:
        sigmaY = sigmaX
        sigmaZ = 0
        experiment = f"x{sigmaX:03}y{sigmaY:03}z{sigmaZ:03}t{sigmaT:03}"
        # Read periodic profile input
        data = np.loadtxt(f"{pathPeriodic}input/prof.inp.000")
        zt   = data[:,0]
        thl0 = data[:,1]
        qt0  = data[:,2]
        u0   = data[:,3]
        v0   = data[:,4]
        e120 = data[:,5]
        ktot = len(zt)

        # Read crossyz data
        nc = Dataset(f"{pathPeriodic}output/crossyz/thlyz.000.nc",'r')
        thlyz = nc.variables['thlyz'][:,:,:]
        yt = nc.variables['yt'][:]
        jtot = len(yt)
        time = nc.variables['time'][:]
        time = time-2*time[0]+time[1]
        nt = len(time)
        nc.close()
        nc = Dataset(f"{pathPeriodic}output/crossyz/qtyz.000.nc",'r')
        qtyz = nc.variables['qtyz'][:,:,:]
        nc.close()
        nc = Dataset(f"{pathPeriodic}output/crossyz/uyz.000.nc",'r')
        uyz = nc.variables['uyz'][:,:,:]
        nc.close()
        nc = Dataset(f"{pathPeriodic}output/crossyz/vyz.000.nc",'r')
        vyz = nc.variables['vyz'][:,:,:]
        ym = nc.variables['ym'][:]
        ym = np.append(ym,ym[-1]+ym[-1]-ym[-2])
        nc.close()
        nc = Dataset(f"{pathPeriodic}output/crossyz/wyz.000.nc",'r')
        wyz = nc.variables['wyz'][:,:,:]
        zm = nc.variables['zm'][:]
        zm = np.append(zm,zm[-1]+(zm[-1]-zm[-2]))
        nc.close()

        # Read crossxz data
        nc = Dataset(f"{pathPeriodic}output/crossxz/thlxz.000.nc",'r')
        thlxz = nc.variables['thlxz'][:,:,:]
        xt = nc.variables['xt'][:]
        itot = len(xt)
        nc.close()
        nc = Dataset(f"{pathPeriodic}output/crossxz/qtxz.000.nc",'r')
        qtxz = nc.variables['qtxz'][:,:,:]
        nc.close()
        nc = Dataset(f"{pathPeriodic}output/crossxz/uxz.000.nc",'r')
        uxz = nc.variables['uxz'][:,:,:]
        xm = nc.variables['xm'][:]
        xm = np.append(xm,xm[-1]+(xm[-1]-xm[-2]))
        nc.close()
        nc = Dataset(f"{pathPeriodic}output/crossxz/vxz.000.nc",'r')
        vxz = nc.variables['vxz'][:,:,:]
        nc.close()
        nc = Dataset(f"{pathPeriodic}output/crossxz/wxz.000.nc",'r')
        wxz = nc.variables['wxz'][:,:,:]
        nc.close()

        # Read crossxy data
        nc = Dataset(f"{pathPeriodic}output/crossxy/{ktot:04}/thlxy.{ktot:04}.000.nc",'r')
        thlxy = nc.variables['thlxy'][:,:,:]
        nc.close()
        nc = Dataset(f"{pathPeriodic}output/crossxy/{ktot:04}/qtxy.{ktot:04}.000.nc",'r')
        qtxy = nc.variables['qtxy'][:,:,:]
        nc.close()
        nc = Dataset(f"{pathPeriodic}output/crossxy/{ktot:04}/uxy.{ktot:04}.000.nc",'r')
        uxy = nc.variables['uxy'][:,:,:]
        nc.close()
        nc = Dataset(f"{pathPeriodic}output/crossxy/{ktot:04}/vxy.{ktot:04}.000.nc",'r')
        vxy = nc.variables['vxy'][:,:,:]
        nc.close()
        wxy = np.zeros((nt,jtot,itot))

        # Downscale spatiotemporal data using a gaussian smoother
        thlxz = sp.gaussian_filter(thlxz,[sigmaT,sigmaZ,sigmaX],order=0,mode=['nearest','nearest','wrap'])
        qtxz = sp.gaussian_filter(qtxz,[sigmaT,sigmaZ,sigmaX],order=0,mode=['nearest','nearest','wrap'])
        uxz = sp.gaussian_filter(uxz,[sigmaT,sigmaZ,sigmaX],order=0,mode=['nearest','nearest','wrap'])
        vxz = sp.gaussian_filter(vxz,[sigmaT,sigmaZ,sigmaX],order=0,mode=['nearest','nearest','wrap'])
        wxz = sp.gaussian_filter(wxz,[sigmaT,sigmaZ,sigmaX],order=0,mode=['nearest','nearest','wrap'])

        thlyz = sp.gaussian_filter(thlyz,[sigmaT,sigmaZ,sigmaY],order=0,mode=['nearest','nearest','wrap'])
        qtyz = sp.gaussian_filter(qtyz,[sigmaT,sigmaZ,sigmaY],order=0,mode=['nearest','nearest','wrap'])
        uyz = sp.gaussian_filter(uyz,[sigmaT,sigmaZ,sigmaY],order=0,mode=['nearest','nearest','wrap'])
        vyz = sp.gaussian_filter(vyz,[sigmaT,sigmaZ,sigmaY],order=0,mode=['nearest','nearest','wrap'])
        wyz = sp.gaussian_filter(wyz,[sigmaT,sigmaZ,sigmaY],order=0,mode=['nearest','nearest','wrap'])

        thlxy = sp.gaussian_filter(thlxy,[sigmaT,sigmaY,sigmaX],order=0,mode=['nearest','wrap','wrap'])
        qtxy = sp.gaussian_filter(qtxy,[sigmaT,sigmaY,sigmaX],order=0,mode=['nearest','wrap','wrap'])
        uxy = sp.gaussian_filter(uxy,[sigmaT,sigmaY,sigmaX],order=0,mode=['nearest','wrap','wrap'])
        vxy = sp.gaussian_filter(vxy,[sigmaT,sigmaY,sigmaX],order=0,mode=['nearest','wrap','wrap'])

        # Add extra half levels
        vyz = np.append(vyz,vyz[:,:,0,None],axis=2)
        uxz = np.append(uxz,uxz[:,:,0,None],axis=2)
        uxy = np.append(uxy,uxy[:,:,0,None],axis=2)
        vxy = np.append(vxy,vxy[:,0,None,:],axis=1)
        wxz = np.append(wxz,np.zeros((nt,1,itot)),axis=1)
        wyz = np.append(wyz,np.zeros((nt,1,jtot)),axis=1)
        wyz[:,0,:] = 0.
        wxz[:,0,:] = 0.

        # Create openboundaries.inp.iexnpr.nc
        openBC=DALES.openBoundariesInp(itot,jtot,ktot,iexpnr=experiment, path=pathWrite)
        # Write dimensions to file
        openBC.writeDimensions(xt,yt,zt,xm,ym,zm)
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
            openBC.writeWest(thlyz[it,:,:],qtyz[it,:,:],uyz[it,:,:],vyz[it,:,:],wyz[it,:,:],e120[:,None]*np.ones((ktot,jtot)),it+1)
            openBC.writeEast(thlyz[it,:,:],qtyz[it,:,:],uyz[it,:,:],vyz[it,:,:],wyz[it,:,:],e120[:,None]*np.ones((ktot,jtot)),it+1)
            openBC.writeSouth(thlxz[it,:,:],qtxz[it,:,:],uxz[it,:,:],vxz[it,:,:],wxz[it,:,:],e120[:,None]*np.ones((ktot,itot)),it+1)
            openBC.writeNorth(thlxz[it,:,:],qtxz[it,:,:],uxz[it,:,:],vxz[it,:,:],wxz[it,:,:],e120[:,None]*np.ones((ktot,itot)),it+1)
            openBC.writeTop(thlxy[it,:,:],qtxy[it,:,:],uxy[it,:,:],vxy[it,:,:],wxy[it,:,:],np.zeros((jtot,itot)),it+1)
        # Close file
        openBC.close()
        del openBC
