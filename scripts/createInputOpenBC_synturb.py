import numpy as np
from netCDF4 import Dataset
import numpy.ma as ma

# Input
pathInput = '../input/boundary_input/'
sigmaX_array = np.array([2])
sigmaT_array = np.array([0])
dxint  = 15360
dyint  = 3840
# Constants
boundarynames = ['west','east','south','north','top']
dim1 = ['zpatch','zpatch','zpatch','zpatch','ypatch']
dim2 = ['ypatch','ypatch','xpatch','xpatch','xpatch']
# Open source data
nc_or = Dataset(f"{pathInput}openboundaries.inp.x000y000z000t000.nc",'r')
xt = nc_or.variables['xt'][:]
yt = nc_or.variables['yt'][:]
zt = nc_or.variables['zt'][:]
time = nc_or.variables['time'][:]
nxpatch = int((xt[1]-xt[0])/dxint*len(xt))
nypatch = int((yt[1]-yt[0])/dyint*len(yt))
nzpatch = len(zt)
xpatch  = np.arange(nxpatch)*dxint+dxint/2
ypatch  = np.arange(nypatch)*dyint+dyint/2
zpatch  = zt

def calc_cov_side(uturb,vturb,wturb,thlturb,qtturb):
    u2 = np.zeros(np.shape(uturb)[:2])
    v2 = np.zeros(np.shape(uturb)[:2])
    w2 = np.zeros(np.shape(uturb)[:2])
    uv = np.zeros(np.shape(uturb)[:2])
    uw = np.zeros(np.shape(uturb)[:2])
    vw = np.zeros(np.shape(uturb)[:2])
    thl2 = np.zeros(np.shape(uturb)[:2])
    wthl = np.zeros(np.shape(uturb)[:2])
    qt2 = np.zeros(np.shape(uturb)[:2])
    wqt = np.zeros(np.shape(uturb)[:2])
    for it in range(np.shape(uturb)[0]):
        for iz in range(np.shape(uturb)[1]):
            covar = np.cov(np.array([ uturb[it,iz,:] , vturb[it,iz,:], wturb[it,iz,:] ]))
            u2[it,iz] = covar[0,0]
            v2[it,iz] = covar[1,1]
            w2[it,iz] = covar[2,2]
            uv[it,iz] = covar[0,1]
            uw[it,iz] = covar[0,2]
            vw[it,iz] = covar[1,2]
            covar = np.cov(np.array([ wturb[it,iz,:], thlturb[it,iz,:] ]))
            thl2[it,iz] = covar[1,1]
            wthl[it,iz] = covar[0,1]
            covar = np.cov(np.array([ wturb[it,iz,:], qtturb[it,iz,:] ]))
            qt2[it,iz] = covar[1,1]
            wqt[it,iz] = covar[0,1] 
    return u2,v2,w2,uv,uw,vw,thl2,wthl,qt2,wqt

def calc_cov_top(uturb,vturb,wturb,thlturb,qtturb):
    u2 = np.zeros(np.shape(uturb)[0])
    v2 = np.zeros(np.shape(uturb)[0])
    w2 = np.zeros(np.shape(uturb)[0])
    uv = np.zeros(np.shape(uturb)[0])
    uw = np.zeros(np.shape(uturb)[0])
    vw = np.zeros(np.shape(uturb)[0])
    thl2 = np.zeros(np.shape(uturb)[0])
    wthl = np.zeros(np.shape(uturb)[0])
    qt2 = np.zeros(np.shape(uturb)[0])
    wqt = np.zeros(np.shape(uturb)[0])
    for it in range(np.shape(uturb)[0]):
        covar = np.cov(np.array([ uturb[it,:,:].flatten() , vturb[it,:,:].flatten(), wturb[it,:,:].flatten() ]))
        u2[it] = covar[0,0]
        v2[it] = covar[1,1]
        w2[it] = covar[2,2]
        uv[it] = covar[0,1]
        uw[it] = covar[0,2]
        vw[it] = covar[1,2]
        covar = np.cov(np.array([ wturb[it,:,:].flatten(), thlturb[it,:,:].flatten() ]))
        thl2[it] = covar[1,1]
        wthl[it] = covar[0,1]
        covar = np.cov(np.array([ wturb[it,:,:].flatten(), qtturb[it,:,:].flatten() ]))
        qt2[it] = covar[1,1]
        wqt[it] = covar[0,1] 
    return u2,v2,w2,uv,uw,vw,thl2,wthl,qt2,wqt

for sigmaX in sigmaX_array:
    for sigmaT in sigmaT_array:
        if(sigmaX==0 and sigmaT==0): continue
        experiment=f"x{sigmaX:03}y{sigmaX:03}z000t{sigmaT:03}"
        print(f"Creating synthetic turbulence input for {experiment}")
        nc = Dataset(f"{pathInput}openboundaries.inp.{experiment}.nc",'r+')
        variables = nc.variables.keys()
        # Add variables to netcdf file
        if(not('zpatch' in variables)):
            nc.createDimension('zpatch',nzpatch)
            nczpatch = nc.createVariable('zpatch','f4',('zpatch'))
            nczpatch[:] = zpatch
        if(not('ypatch' in variables)):
            nc.createDimension('ypatch',nypatch)
            ncypatch = nc.createVariable('ypatch','f4',('ypatch'))
            ncypatch[:] = ypatch
        if(not('xpatch' in variables)):
            nc.createDimension('xpatch',nxpatch)
            ncxpatch = nc.createVariable('xpatch','f4',('xpatch'))
            ncxpatch[:] = xpatch
        for ib in range(len(boundarynames)):
            boundaryname = boundarynames[ib]
            if(not(f"u2{boundaryname}" in variables)):
                u2   = nc.createVariable(f"u2{boundaryname}",'f4',('time',dim1[ib],dim2[ib]))
            else:
                u2   = nc.variables[f"u2{boundaryname}"]
            if(not(f"v2{boundaryname}" in variables)):
                v2   = nc.createVariable(f"v2{boundaryname}",'f4',('time',dim1[ib],dim2[ib]))
            else:
                v2   = nc.variables[f"v2{boundaryname}"]
            if(not(f"w2{boundaryname}" in variables)):
                w2   = nc.createVariable(f"w2{boundaryname}",'f4',('time',dim1[ib],dim2[ib]))
            else:
                w2   = nc.variables[f"w2{boundaryname}"]
            if(not(f"uv{boundaryname}" in variables)):
                uv   = nc.createVariable(f"uv{boundaryname}",'f4',('time',dim1[ib],dim2[ib]))
            else:
                uv   = nc.variables[f"uv{boundaryname}"]
            if(not(f"uw{boundaryname}" in variables)):
                uw   = nc.createVariable(f"uw{boundaryname}",'f4',('time',dim1[ib],dim2[ib]))
            else:
                uw   = nc.variables[f"uw{boundaryname}"]
            if(not(f"vw{boundaryname}" in variables)):
                vw   = nc.createVariable(f"vw{boundaryname}",'f4',('time',dim1[ib],dim2[ib]))
            else:
                vw   = nc.variables[f"vw{boundaryname}"]
            if(not(f"thl2{boundaryname}" in variables)):
                thl2 = nc.createVariable(f"thl2{boundaryname}",'f4',('time',dim1[ib],dim2[ib]))
            else:
                thl2 = nc.variables[f"thl2{boundaryname}"]
            if(not(f"wthl{boundaryname}" in variables)):
                wthl = nc.createVariable(f"wthl{boundaryname}",'f4',('time',dim1[ib],dim2[ib]))
            else:
                wthl = nc.variables[f"wthl{boundaryname}"]
            if(not(f"qt2{boundaryname}" in variables)):
                qt2  = nc.createVariable(f"qt2{boundaryname}",'f4',('time',dim1[ib],dim2[ib]))
            else:
                qt2  = nc.variables[f"qt2{boundaryname}"]
            if(not(f"wqt{boundaryname}" in variables)):
                wqt  = nc.createVariable(f"wqt{boundaryname}",'f4',('time',dim1[ib],dim2[ib]))
            else:
                wqt  = nc.variables[f"wqt{boundaryname}"]
            # Read original data
            u  = nc_or.variables[f"u{boundaryname}"][:,:,:]
            v  = nc_or.variables[f"v{boundaryname}"][:,:,:]
            w  = nc_or.variables[f"w{boundaryname}"][:,:,:]
            thl = nc_or.variables[f"thl{boundaryname}"][:,:,:]
            qt  = nc_or.variables[f"qt{boundaryname}"][:,:,:]
            # Read filtered data
            umean   = nc.variables[f"u{boundaryname}"][:,:,:]
            vmean   = nc.variables[f"v{boundaryname}"][:,:,:]
            wmean   = nc.variables[f"w{boundaryname}"][:,:,:]
            thlmean = nc.variables[f"thl{boundaryname}"][:,:,:]
            qtmean  = nc.variables[f"qt{boundaryname}"][:,:,:]
            # Get all variables at grid centres
            if(ib==0 or ib==1):
                v = (v[:,:,1:]+v[:,:,:-1])/2
                w = (w[:,1:,:]+w[:,:-1,:])/2
                vmean = (vmean[:,:,1:]+vmean[:,:,:-1])/2
                wmean = (wmean[:,1:,:]+wmean[:,:-1,:])/2
            elif(ib==2 or ib==3):
                u = (u[:,:,1:]+u[:,:,:-1])/2
                w = (w[:,1:,:]+w[:,:-1,:])/2
                umean = (umean[:,:,1:]+umean[:,:,:-1])/2
                wmean = (wmean[:,1:,:]+wmean[:,:-1,:])/2
            else:
                u = (u[:,:,1:]+u[:,:,:-1])/2
                v = (v[:,1:,:]+v[:,:-1,:])/2
                umean = (umean[:,:,1:]+umean[:,:,:-1])/2
                vmean = (vmean[:,1:,:]+vmean[:,:-1,:])/2
            # Obtain perturbation fields
            uturb   = ma.getdata(umean-u)
            vturb   = ma.getdata(vmean-v)
            wturb   = ma.getdata(wmean-w)
            thlturb = ma.getdata(thlmean-thl)
            qtturb  = ma.getdata(qtmean-qt)
            if(ib!=4):
                u2data,v2data,w2data,uvdata,uwdata,vwdata,thl2data,wthldata,qt2data,wqtdata = \
                    calc_cov_side(uturb,vturb,wturb,thlturb,qtturb)
                u2[:] = u2data
                v2[:] = v2data
                w2[:] = w2data
                uv[:] = uvdata
                uw[:] = uwdata
                vw[:] = vwdata
                thl2[:] = thl2data
                wthl[:] = wthldata
                qt2[:] = qt2data
                wqt[:] = wqtdata
            else:
                u2data,v2data,w2data,uvdata,uwdata,vwdata,thl2data,wthldata,qt2data,wqtdata = \
                    calc_cov_top(uturb,vturb,wturb,thlturb,qtturb)
                u2[:] = u2data
                v2[:] = v2data
                w2[:] = w2data
                uv[:] = uvdata
                uw[:] = uwdata
                vw[:] = vwdata
                thl2[:] = thl2data
                wthl[:] = wthldata
                qt2[:] = qt2data
                wqt[:] = wqtdata
        nc.close()
nc_or.close()
