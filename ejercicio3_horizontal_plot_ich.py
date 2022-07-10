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
archive = 'ich_bio.nc'
# --------------------------------------------------------------------------+

particle_size = 1

# Mapping selection zone with lon and lat of the globe

# Atlantic coastline
minlon = -12
maxlon = 0
minlat = 36.5
maxlat = 46.5
initime = datetime.datetime(2006,1,1)



# --------------------------------------------------------------------------+

nc = netCDF4.Dataset(archive)
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
z = nc.variables['depth'][:]
time = nc.variables['time'][:]


# --------------------------------------------------------------------------+

fintime=[]

for t in time:
    fintime.append(initime+datetime.timedelta(t/86400.))
    
dates_list = [datetime.datetime(2006, 11, 15, 0, 0), datetime.datetime(2006, 11, 20, 0, 0), datetime.datetime(2006, 11, 25, 0, 0), datetime.datetime(2006, 12, 1, 0, 0)]

index_list = []

for i in range(len(dates_list)):
    for j in range(len(fintime)):
        if (fintime[j] == dates_list[i]) == True:
            index_list.append(j) 

# --------------------------------------------------------------------------+

fig = plt.figure(figsize = (10,10))
map = Basemap(llcrnrlon = minlon, llcrnrlat = minlat, urcrnrlon = maxlon, urcrnrlat = maxlat, resolution = 'f',  projection = 'merc', lat_0 = minlat, lon_0 = minlon)

Lon, Lat = map(lon, lat)




plt.plot(Lon, Lat, color = 'r', marker = '.', linewidth = 0, markersize = particle_size)
#plt.plot(Lon_m[:,0], Lat_m[:,0], color = 'b', marker = '.', linewidth = 0, markersize = particle_size)
plt.plot(Lon[-1,:], Lat[-1,:], color = 'g', marker = '.', linewidth = 0, markersize = particle_size)

#plt.text(np.mean(Lon[:,0]),np.mean(Lat[:,0]),'Lanzamento',fontweight='bold')
	
# --------------------------------------------------------------------------+ 
   
map.fillcontinents(color = 'grey')
parallels = np.arange(minlat, maxlat, 2)
meridians = np.arange(minlon, maxlon, 2)

map.drawparallels(parallels, labels = [1,0,0,0], fontsize = 10)
map.drawmeridians(meridians, labels = [0,0,0,1], fontsize = 10)


plt.show()

fig.savefig('ich_total.png')

# --------------------------------------------------------------------------+ 

fig_day_1 = plt.figure(figsize = (10,10))

plt.plot(Lon[index_list[0],:], Lat[index_list[0],:], color = 'r', marker = '.', linewidth = 0, markersize = particle_size)
plt.title(dates_list[0].strftime("%d-%m-%Y, %H:%M:%S"))

map.fillcontinents(color = 'grey')
parallels = np.arange(minlat, maxlat, 2)
meridians = np.arange(minlon, maxlon, 2)

map.drawparallels(parallels, labels = [1,0,0,0], fontsize = 10)
map.drawmeridians(meridians, labels = [0,0,0,1], fontsize = 10)

plt.show()

fig_day_1.savefig('ich_total_15.png')

# --------------------------------------------------------------------------+ 

fig_day_2 = plt.figure(figsize = (10,10))

plt.plot(Lon[index_list[1],:], Lat[index_list[1],:], color = 'r', marker = '.', linewidth = 0, markersize = particle_size)
plt.title(dates_list[1].strftime("%d-%m-%Y, %H:%M:%S"))

map.fillcontinents(color = 'grey')
parallels = np.arange(minlat, maxlat, 2)
meridians = np.arange(minlon, maxlon, 2)

map.drawparallels(parallels, labels = [1,0,0,0], fontsize = 10)
map.drawmeridians(meridians, labels = [0,0,0,1], fontsize = 10)

plt.show()

fig_day_2.savefig('ich_total_20.png')

# --------------------------------------------------------------------------+ 

fig_day_3 = plt.figure(figsize = (10,10))

plt.plot(Lon[index_list[2],:], Lat[index_list[2],:], color = 'r', marker = '.', linewidth = 0, markersize = particle_size)
plt.title(dates_list[2].strftime("%d-%m-%Y, %H:%M:%S"))

map.fillcontinents(color = 'grey')
parallels = np.arange(minlat, maxlat, 2)
meridians = np.arange(minlon, maxlon, 2)

map.drawparallels(parallels, labels = [1,0,0,0], fontsize = 10)
map.drawmeridians(meridians, labels = [0,0,0,1], fontsize = 10)

plt.show()

fig_day_3.savefig('ich_total_25.png')

# --------------------------------------------------------------------------+ 

fig_day_4 = plt.figure(figsize = (10,10))

plt.plot(Lon[index_list[3],:], Lat[index_list[3],:], color = 'r', marker = '.', linewidth = 0, markersize = particle_size)
plt.title(dates_list[3].strftime("%d-%m-%Y, %H:%M:%S"))

map.fillcontinents(color = 'grey')
parallels = np.arange(minlat, maxlat, 2)
meridians = np.arange(minlon, maxlon, 2)

map.drawparallels(parallels, labels = [1,0,0,0], fontsize = 10)
map.drawmeridians(meridians, labels = [0,0,0,1], fontsize = 10)

plt.show()

fig_day_4.savefig('ich_total_30.png')
