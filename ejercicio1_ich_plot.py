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

# --------------------------------------------------------------------------+


archive = 'output.nc'



particle_size = 1

minlon = -10
maxlon = -8.5
minlat = 39
maxlat = 41
initime = datetime.datetime(2005,1,1)


# --------------------------------------------------------------------------+

nc = netCDF4.Dataset(archive)
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
z = nc.variables['depth'][:]
time = nc.variables['time'][:]
mortality = nc.variables['mortality'][:]



# --------------------------------------------------------------------------+

fig = plt.figure(figsize = (10,10))
map = Basemap(llcrnrlon = minlon, llcrnrlat = minlat, urcrnrlon = maxlon, urcrnrlat = maxlat, resolution = 'f',  projection = 'merc', lat_0 = minlat, lon_0 = minlon)

Lon,Lat=map(lon,lat)

Lon_alive = Lon.copy()
Lat_alive = Lat.copy()

for i in range(1, len(Lon_alive[:,0])):
    Lon_alive[i,:] = ma.masked_where(mortality[i,:] == 4, Lon_alive[i,:])
    Lat_alive[i,:] = ma.masked_where(mortality[i,:] == 4, Lat_alive[i,:])

plt.plot(Lon, Lat, color = 'r', marker = '.', linewidth = 0, markersize = particle_size)
plt.plot(Lon[0,:], Lat[0,:], color = 'b', marker = '.', linewidth = 0, markersize = particle_size)
plt.plot(Lon[-1,:], Lat[-1,:], color = 'g', marker = '.', linewidth = 0, markersize = particle_size)

plt.text(np.mean(Lon[0,:]),np.mean(Lat[0,:]),'Lanzamento',fontweight='bold')
	
# --------------------------------------------------------------------------+ 
   
map.fillcontinents(color = 'grey')
parallels = np.arange(minlat, maxlat, 0.5)
meridians = np.arange(minlon, maxlon, 0.5)

map.drawparallels(parallels, labels = [1,0,0,0], fontsize = 10)
map.drawmeridians(meridians, labels = [0,0,0,1], fontsize = 10)

namefig = 'Lanzamientos.jpg'


plt.show()
fig.savefig('horizontal.png')

# gc

fig_mean = plt.figure(figsize = (10,10))
map = Basemap(llcrnrlon = minlon, llcrnrlat = minlat, urcrnrlon = maxlon, urcrnrlat = maxlat, resolution = 'f',  projection = 'merc', lat_0 = minlat, lon_0 = minlon)

fintime=[]
for t in time:
    fintime.append(initime+datetime.timedelta(t/86400.))

plt.plot(Lon, Lat, color = 'r', marker = '.', linewidth = 0, markersize = particle_size)
for i in range(len(fintime)):
    plt.plot(np.mean(Lon_alive[i,:]),np.mean(Lat_alive[i,:]), color = 'b', marker = '.', linewidth = 0, markersize = particle_size)

plt.text(np.mean(Lon[0,:]),np.mean(Lat[0,:]),'Lanzamento',fontweight='bold')

map.fillcontinents(color = 'grey')
parallels = np.arange(minlat, maxlat, 0.5)
meridians = np.arange(minlon, maxlon, 0.5)

map.drawparallels(parallels, labels = [1,0,0,0], fontsize = 10)
map.drawmeridians(meridians, labels = [0,0,0,1], fontsize = 10)

plt.show()

fig_mean.savefig('gc.png')
