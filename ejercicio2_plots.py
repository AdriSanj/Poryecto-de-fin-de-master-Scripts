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
phase = 2

figname = 'At_Ah'
archive1 = figname + '_ich.nc'
archive2 = figname + '_od.nc'



particle_size = 1

#minlon = -10.5
#maxlon = -8
#minlat = 39
#maxlat = 42
minlon = -12.5
maxlon = -11.5
minlat = 42.5
maxlat = 43.5
initime = datetime.datetime(2009,1,1)



nc1 = netCDF4.Dataset(archive1)
lon1 = nc1.variables['lon'][:]
lat1 = nc1.variables['lat'][:]
z1 = nc1.variables['depth'][:]
time1 = nc1.variables['time'][:]
mortality1 = nc1.variables['mortality'][:]


nc2 = netCDF4.Dataset(archive2)
lon2 = nc2.variables['lon'][:]
lat2 = nc2.variables['lat'][:]
z2 = nc2.variables['z'][:]
time2 = nc2.variables['time'][:]


# --------------------------------------------------------------------------+
initime_ich=datetime.datetime(1900,1,1)
fintime_ich=[]
for t in time1:
    fintime_ich.append(initime_ich+datetime.timedelta(t/86400.))
    
initime_od=datetime.datetime(1970,1,1)
fintime_od=[]
for t in time2:
    fintime_od.append(initime_od+datetime.timedelta(t/86400.))
    
map = Basemap(llcrnrlon = minlon, llcrnrlat = minlat, urcrnrlon = maxlon, urcrnrlat = maxlat, resolution = 'f',  projection = 'merc', lat_0 = minlat, lon_0 = minlon)

Lon1,Lat1=map(lon1,lat1)

Lon_alive1 = Lon1.copy()
Lat_alive1 = Lat1.copy()

Lon2,Lat2=map(lon2,lat2)

Lon_all2 = Lon2.copy()
Lat_all2 = Lat2.copy()

for i in range(1, len(Lon_alive1[:,0])):
    Lon_alive1[i,:] = ma.masked_where(mortality1[i,:] == 4, Lon_alive1[i,:])
    Lat_alive1[i,:] = ma.masked_where(mortality1[i,:] == 4, Lat_alive1[i,:])
    
for i in range(1, len(Lon_all2[0,:])):
    Lon_all2[:,i] = np.ma.filled(Lon_all2[:,i].astype(float), Lon_all2[:,i-1])
    Lat_all2[:,i] = np.ma.filled(Lat_all2[:,i].astype(float), Lat_all2[:,i-1])
    
    
z1_m = z1.copy()

for i in range(1, len(fintime_ich)):
    z1_m[i,:] = ma.masked_where(mortality1[i,:] == 4, z1_m[i,:])


fig1 = plt.figure(figsize = (10,10))
plt.plot(Lon1, Lat1, color = 'r', marker = '.', linewidth = 0, markersize = particle_size)
plt.plot(Lon1[0,:], Lat1[0,:], color = 'b', marker = '.', linewidth = 0, markersize = particle_size)
plt.plot(Lon1[-1,:], Lat1[-1,:], color = 'g', marker = '.', linewidth = 0, markersize = particle_size)

plt.text(np.mean(Lon1[0,:]),np.mean(Lat1[0,:]),'Lanzamento',fontweight='bold')
	
# --------------------------------------------------------------------------+ 
   
map.fillcontinents(color = 'grey')
parallels = np.arange(minlat, maxlat, 0.5)
meridians = np.arange(minlon, maxlon, 0.5)

map.drawparallels(parallels, labels = [1,0,0,0], fontsize = 10)
map.drawmeridians(meridians, labels = [0,0,0,1], fontsize = 10)



plt.show()
fig1.savefig(figname+'_ich.png')



fig2 = plt.figure(figsize = (10,10))
plt.plot(Lon_all2, Lat_all2, color = 'r', marker = '.', linewidth = 0, markersize = particle_size)
plt.plot(Lon_all2[:,0], Lat_all2[:,0], color = 'b', marker = '.', linewidth = 0, markersize = particle_size)
plt.plot(Lon_all2[:,-1], Lat_all2[:,-1], color = 'g', marker = '.', linewidth = 0, markersize = particle_size)

plt.text(np.mean(Lon_all2[:,0]),np.mean(Lat_all2[:,0]),'Lanzamento',fontweight='bold')
	
# --------------------------------------------------------------------------+ 
   
map.fillcontinents(color = 'grey')
parallels = np.arange(minlat, maxlat, 0.5)
meridians = np.arange(minlon, maxlon, 0.5)

map.drawparallels(parallels, labels = [1,0,0,0], fontsize = 10)
map.drawmeridians(meridians, labels = [0,0,0,1], fontsize = 10)

plt.show()

fig2.savefig(figname+'_od.png')

# gc

if phase > 1:
    fig3 = plt.figure(figsize = (10,10))
    
    
    for i in range(len(fintime_ich)):
        plt.plot(np.mean(Lon_alive1[i,:]),np.mean(Lat_alive1[i,:]), color = 'r', marker = '.', linewidth = 0, markersize = particle_size)
        plt.plot(np.mean(Lon2[:,i]),np.mean(Lat2[:,i]), color = 'b', marker = '.', linewidth = 0, markersize = particle_size)
        
    plt.text(np.mean(Lon1[0,:]),np.mean(Lat1[0,:]),'Lanzamento',fontweight='bold')

    map.fillcontinents(color = 'grey')
    parallels = np.arange(minlat, maxlat, 0.5)
    meridians = np.arange(minlon, maxlon, 0.5)

    map.drawparallels(parallels, labels = [1,0,0,0], fontsize = 10)
    map.drawmeridians(meridians, labels = [0,0,0,1], fontsize = 10)

    plt.show()

    fig3.savefig(figname+'_gc.png')







if phase > 2:   
    medias_ICH = []
    medias_OD_alive = []
    for i in range(len(fintime_od)):
        medias_OD_alive.append(np.mean(z2[:,i]))
    for j in range(len(fintime_ich)):
        medias_ICH.append(np.mean(z1_m[j,:]))
    
    fig4 = plt.figure()

    plt.plot(fintime_od,medias_OD_alive,'b-')
    plt.plot(fintime_ich, medias_ICH,'r-')
    plt.title('Profundidad vs Tiempo')
    plt.xlabel('Tiempo')
    plt.ylabel('Profundidad (m)')
    plt.grid()
    plt.show()
    fig4.savefig(figname+'_mean_depths.png')
