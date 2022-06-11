
import sys
sys.path.append("/home/xuewei/3rdParty/anaconda3/pkgs/lib/python3.6/site-packages")
import cmaps
import matplotlib
matplotlib.use('Agg')
import numpy as np
from mpi4py import MPI
import netCDF4 as nc
import f90nml
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

comm = MPI.COMM_WORLD
pid  = comm.Get_rank()


namelistpath = sys.argv[1]
nml      = f90nml.read(namelistpath)
nprocs   = int(nml['nml']['nprocs'])
ix       = int(nml['nml']['ix'])
iy       = int(nml['nml']['iy'])
iz       = int(nml['nml']['iz'])
var      = nml['nml']['var'].split(',')
member   = int(nml['nml']['member'])
casename = nml['nml']['casename']
fpath    = nml['nml']['fpath']


fid = nc.Dataset(fpath)
tinterval= int(1024/member)

num_lat = fid.dimensions['lat'].size
num_lon = fid.dimensions['lon'].size
num_lev = fid.dimensions['lev'].size
num_ilat= fid.dimensions['ilat'].size
num_ilon= fid.dimensions['ilon'].size


lon = np.zeros((num_lat, num_lon))
lat = np.zeros((num_lat, num_lon))
ilon= np.zeros((num_lat, num_ilon))
ilat= np.zeros((num_ilat, num_lon))

lon[0,:] = fid.variables['lon'][:]
for i in range(1, num_lat):
    lon[i,:] = lon[0,:]
lat[:,0] = fid.variables['lat'][:]
for i in range(1, num_lon):
    lat[:,i] = lat[:,0]
ilon[0,:] = fid.variables['ilon'][:]
for i in range(1, num_lat):
    ilon[i,:] = ilon[0,:]
ilat[:,0] = fid.variables['ilat'][:]
for i in range(1, num_lon):
    ilat[:,i] = ilat[:,0]


print(lon)
sec = num_lon // nprocs
res = np.mod(num_lon, nprocs)
print ('a',sec,res)
if ( pid < res):
    start = int(pid * (sec+1))
    end   = int(start + sec)
else:
    start = int(pid * sec + res)
    end   = int(start + sec - 1)

for v in var:
    if v == 'phs':
        ref_value= fid.variables[v][0:1023:tinterval, iy-1, ix-1]
        value = fid.variables[v][0:1023:tinterval, :, start:end+1]
        mean_ref_value = np.mean(ref_value, 0)
        mean_value= np.mean(value, 0)
        for t in range(ref_value.shape[0]):
            ref_value[t] = ref_value[t] - mean_ref_value
            value[t,:,:] = value[t,:,:] - mean_value
    else:
        ref_value= fid.variables[v][0:1023:tinterval, iz-1, iy-1, ix-1]
        value = fid.variables[v][0:1023:tinterval, iz-1, :, start:end+1]
        mean_ref_value = np.mean(ref_value, 0)
        mean_value= np.mean(value, 0)
        for t in range(ref_value.shape[0]):
            ref_value[t] = ref_value[t] - mean_ref_value
            value[t,:,:] = value[t,:,:] - mean_value

    if v == 'v':
        corr = np.zeros((num_ilat, end-start+1))
        num_lat_tmp = num_ilat
    else:
        corr = np.zeros((num_lat, end-start+1))
        num_lat_tmp = num_lat

    for i in range(end-start+1):
        for j in range(num_lat_tmp):
            if (np.var(value[:,j,i])==0):
                corr[j, i] = None
                continue
            tmp = np.corrcoef(ref_value, value[:,j,i])
            corr[j, i] = tmp[0,1]
    
    if pid != 0:
        comm.send(corr, dest=0)
    else:
        if v == 'v':
            corr_all = np.zeros((num_ilat, num_lon))
        else:
            corr_all = np.zeros((num_lat, num_lon))
        corr_all[:,start:end+1] = corr
        for n in range(1, nprocs):
            if ( n < res):
                rstart = int(n * (sec+1))
                rend   = int(rstart + sec)
            else:
                rstart = int(n * sec + res)
                rend   = int(rstart + sec - 1)
            print (rstart, rend)
            corr_all[:,rstart:rend+1] = comm.recv(source=n)
        
        fig = plt.figure(figsize=(20, 15))
        basemap = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=360,resolution='c')

        if v == 'u':
            pcl = basemap.pcolormesh(ilon, lat, corr_all, shading='gouraud', cmap=cmaps.BlueWhiteOrangeRed, vmin=-1, vmax=1)
            CS = basemap.contour(ilon, lat, corr_all, colors='black')
        elif v == 'v':
            pcl = basemap.pcolormesh(lon[0:num_ilat,:], ilat, corr_all, shading='gouraud', cmap=cmaps.BlueWhiteOrangeRed, vmin=-1, vmax=1)
            CS = basemap.contour(lon[0:num_ilat,:], ilat, corr_all, colors='black')
        else:
            pcl = basemap.pcolormesh(lon, lat, corr_all, shading='gouraud', cmap=cmaps.BlueWhiteOrangeRed, vmin=-1, vmax=1)
            CS = basemap.contour(lon, lat, corr_all, colors='black')



        basemap.scatter(lon[0,ix-1], lat[iy-1,0], s=200, c='y', marker='*')
        basemap.drawparallels(np.arange(-90.,91.,20.), labels=[1,0,0,0], fontsize=20, linewidth=0)
        basemap.drawmeridians(np.arange(0., 361., 40.), labels=[0,0,0,1], fontsize=20, linewidth=0)

        cbar = plt.colorbar(pcl, orientation='horizontal', pad=0.1, aspect=20) 
        cbar.ax.tick_params(labelsize=20) 
        plt.title(v+'-corr-horizontal-'+str(member)+'mem-lev'+str(iz), fontsize=20) 
        plt.savefig(casename+'-'+v+'-corr-horizontal-'+str(member)+'mem-lev'+str(iz)+".png")



