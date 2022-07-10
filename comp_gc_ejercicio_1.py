import numpy as np
import sys
#import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4
import datetime
import glob
import os
import numpy.ma as ma
os.environ['PROJ_LIB'] ='/opt/miniconda/miniconda3/share/proj'
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata


archive_1 = 'ich_ieo.nc'

archive_2 = 'ich_meteo.nc'

archive_3 = 'od_meteo.nc'

archive_4 = 'od_ieo.nc'

particle_size = 1

minlon = -10
maxlon = -8.5
minlat = 39
maxlat = 40.5

# -----------------------------------

nc_1 = netCDF4.Dataset(archive_1)
lon_1 = nc_1.variables['lon'][:]
lat_1 = nc_1.variables['lat'][:]
time = nc_1.variables['time'][:]
mortality_1 = nc_1.variables['mortality'][:]
initime = datetime.datetime(2005,1,1)
fintime=[]
for t in time:
    fintime.append(initime+datetime.timedelta(t/86400.))
    
# -----------------------------------
nc_2 = netCDF4.Dataset(archive_2)
lon_2 = nc_2.variables['lon'][:]
lat_2 = nc_2.variables['lat'][:]
mortality_2 = nc_2.variables['mortality'][:]
# -----------------------------------

nc_3 = netCDF4.Dataset(archive_3)
lon_3 = nc_3.variables['lon'][:]
lat_3 = nc_3.variables['lat'][:]

nc_4 = netCDF4.Dataset(archive_4)
lon_4 = nc_4.variables['lon'][:]
lat_4 = nc_4.variables['lat'][:]

fig = plt.figure(figsize = (10,10))
map = Basemap(llcrnrlon = minlon, llcrnrlat = minlat, urcrnrlon = maxlon, urcrnrlat = maxlat, resolution = 'f',  projection = 'merc', lat_0 = minlat, lon_0 = minlon)

Lon_1,Lat_1=map(lon_1,lat_1)
Lon_2,Lat_2=map(lon_2,lat_2)
Lon_3,Lat_3=map(lon_3,lat_3)
Lon_4,Lat_4=map(lon_4,lat_4)

Lon_alive_1 = Lon_1.copy()
Lat_alive_1 = Lat_1.copy()
Lon_alive_2 = Lon_2.copy()
Lat_alive_2 = Lat_2.copy()

for i in range(1, len(Lon_alive_1[:,0])):
    Lon_alive_1[i,:] = ma.masked_where(mortality_1[i,:] == 4, Lon_alive_1[i,:])
    Lat_alive_1[i,:] = ma.masked_where(mortality_1[i,:] == 4, Lat_alive_1[i,:])
    Lon_alive_2[i,:] = ma.masked_where(mortality_2[i,:] == 4, Lon_alive_2[i,:])
    Lat_alive_2[i,:] = ma.masked_where(mortality_2[i,:] == 4, Lat_alive_2[i,:])


for i in range(len(fintime)):
    plt.plot(np.mean(Lon_alive_1[i,:]),np.mean(Lat_alive_1[i,:]), color = 'b', marker = '.', linewidth = 0, markersize = particle_size)
    plt.plot(np.mean(Lon_alive_2[i,:]),np.mean(Lat_alive_2[i,:]), color = 'r', marker = '.', linewidth = 0, markersize = particle_size)
    plt.plot(np.mean(Lon_3[:,i]),np.mean(Lat_3[:,i]), color = 'g', marker = '.', linewidth = 0, markersize = particle_size)
    plt.plot(np.mean(Lon_4[:,i]),np.mean(Lat_4[:,i]), color = 'y', marker = '.', linewidth = 0, markersize = particle_size)

plt.text(np.mean(Lon_1[0,:]),np.mean(Lat_1[0,:]),'Lanzamento',fontweight='bold')



map.fillcontinents(color = 'grey')
parallels = np.arange(minlat, maxlat, 0.5)
meridians = np.arange(minlon, maxlon, 0.5)

map.drawparallels(parallels, labels = [1,0,0,0], fontsize = 10)
map.drawmeridians(meridians, labels = [0,0,0,1], fontsize = 10)

plt.show()

fig.savefig('gc_comparative.png')

