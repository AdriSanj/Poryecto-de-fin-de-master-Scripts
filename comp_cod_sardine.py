import numpy as np
import sys
#import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4
import datetime
import glob
import os
os.environ['PROJ_LIB'] ='/opt/miniconda/miniconda3/share/proj'
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata

# --------------------------------------------------------------------------+
archive_1 = 'cod_output.nc'
archive_2 = 'sardine_output.nc'
# --------------------------------------------------------------------------+

particle_size = 1

# Mapping selection zone with lon and lat of the globe

# Atlantic coastline
minlon = -10
maxlon = 0
minlat = 41
maxlat = 47

initime = datetime.datetime(1970,1,1)

# --------------------------------------------------------------------------+

nc1 = netCDF4.Dataset(archive_1)
lon1 = nc1.variables['lon'][:]
lat1 = nc1.variables['lat'][:]
z1 = nc1.variables['z'][:]
time1 = nc1.variables['time'][:]
l1 = nc1.variables['length'][:]
s_f1 = nc1.variables['stage_fraction']


# --------------------------------------------------------------------------+

nc2 = netCDF4.Dataset(archive_2)
lon2 = nc2.variables['lon'][:]
lat2 = nc2.variables['lat'][:]
z2 = nc2.variables['z'][:]
time2 = nc2.variables['time'][:]
l2 = nc2.variables['length'][:]
s_f2 = nc2.variables['stage_fraction']


# --------------------------------------------------------------------------+

# --------------------------------------------------------------------------+

fintime=[]

for t in time1:
    fintime.append(initime+datetime.timedelta(t/86400.))
    

# --------------------------------------------------------------------------+


map = Basemap(llcrnrlon = minlon, llcrnrlat = minlat, urcrnrlon = maxlon, urcrnrlat = maxlat, resolution = 'f',  projection = 'merc', lat_0 = minlat, lon_0 = minlon)

Lon1, Lat1 = map(lon1, lat1)
Lon2, Lat2 = map(lon2, lat2)

# Centros de gravedad de las particulas vivas
Lon_m1 = Lon1.copy()
Lat_m1 = Lat1.copy()

Lon_m2 = Lon2.copy()
Lat_m2 = Lat2.copy()

for i in range(1, len(Lon1[0,:])):
    Lon_m1[:,i] = np.ma.filled(Lon_m1[:,i].astype(float), Lon_m1[:,i-1])
    Lat_m1[:,i] = np.ma.filled(Lat_m1[:,i].astype(float), Lat_m1[:,i-1])
    Lon_m2[:,i] = np.ma.filled(Lon_m2[:,i].astype(float), Lon_m2[:,i-1])
    Lat_m2[:,i] = np.ma.filled(Lat_m2[:,i].astype(float), Lat_m2[:,i-1])


medias_OD_1 = []
medias_length_1 = []
medias_OD_2 = []
medias_length_2 = []

for i in range(len(fintime)):
    medias_OD_1.append(np.mean(z1[:,i]))
    medias_length_1.append(np.mean(l1[:,i]))
    medias_OD_2.append(np.mean(z2[:,i]))
    medias_length_2.append(np.mean(l2[:,i]))

fig_1 = plt.figure(figsize = (10,10))

plt.plot(Lon1, Lat1, color = 'r', marker = '.', linewidth = 0, markersize = particle_size)
#plt.plot(Lon_m[:,0], Lat_m[:,0], color = 'b', marker = '.', linewidth = 0, markersize = particle_size)
plt.plot(Lon_m1[:,-1], Lat_m1[:,-1], color = 'g', marker = '.', linewidth = 0, markersize = particle_size)

plt.text(np.mean(Lon1[:,0]),np.mean(Lat1[:,0]),'Lanzamento',fontweight='bold')
	
   
map.fillcontinents(color = 'grey')
parallels = np.arange(minlat, maxlat, 2)
meridians = np.arange(minlon, maxlon, 2)

map.drawparallels(parallels, labels = [1,0,0,0], fontsize = 10)
map.drawmeridians(meridians, labels = [0,0,0,1], fontsize = 10)


plt.show()
fig_1.savefig('horzontal_cod.png')

fig_2 = plt.figure(figsize = (10,10))

plt.plot(Lon2, Lat2, color = 'r', marker = '.', linewidth = 0, markersize = particle_size)
#plt.plot(Lon_m[:,0], Lat_m[:,0], color = 'b', marker = '.', linewidth = 0, markersize = particle_size)
plt.plot(Lon_m2[:,-1], Lat_m2[:,-1], color = 'g', marker = '.', linewidth = 0, markersize = particle_size)

plt.text(np.mean(Lon2[:,0]),np.mean(Lat2[:,0]),'Lanzamento',fontweight='bold')
	
# --------------------------------------------------------------------------+ 
   
map.fillcontinents(color = 'grey')
parallels = np.arange(minlat, maxlat, 2)
meridians = np.arange(minlon, maxlon, 2)

map.drawparallels(parallels, labels = [1,0,0,0], fontsize = 10)
map.drawmeridians(meridians, labels = [0,0,0,1], fontsize = 10)


plt.show()
fig_2.savefig('horizontal_sardine.png')
# --------------------------------------------------------------------------+ 
fig_3 = plt.figure(figsize = (10,10))

for i in range(len(fintime)):
    plt.plot(np.mean(Lon1[:,i]),np.mean(Lat1[:,i]), color = 'b', marker = '.', linewidth = 0, markersize = particle_size)
    plt.plot(np.mean(Lon2[:,i]),np.mean(Lat2[:,i]), color = 'r', marker = '.', linewidth = 0, markersize = particle_size)

plt.text(np.mean(Lon2[:,0]),np.mean(Lat2[:,0]),'Lanzamento',fontweight='bold')

map.fillcontinents(color = 'grey')
parallels = np.arange(minlat, maxlat, 2)
meridians = np.arange(minlon, maxlon, 2)

map.drawparallels(parallels, labels = [1,0,0,0], fontsize = 10)
map.drawmeridians(meridians, labels = [0,0,0,1], fontsize = 10)


plt.show()
fig_3.savefig('gc.png')

# Depth
fig_4 = plt.figure()
plt.plot(fintime, medias_OD_1, 'b')
plt.plot(fintime, medias_OD_2, 'r')
plt.xlabel('Tiempo')
plt.ylabel('Profundidad (m)')
plt.title('Profundidad frente a tiempo')
plt.grid()
plt.show()

fig_4.savefig('mean_depth_cod_sardine.png')

fig_5 = plt.figure()

# Length

plt.plot(fintime, medias_length_1, 'b')
plt.plot(fintime, medias_length_2, 'r')
plt.xlabel('Tiempo')
plt.ylabel('Longitud (mm)')
plt.title('Longitud frente a tiempo')
plt.grid()
plt.show()
fig_5.savefig('mean_length_cod_sardine.png')

