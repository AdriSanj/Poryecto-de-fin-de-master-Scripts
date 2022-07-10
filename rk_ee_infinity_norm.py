import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import datetime
import numpy.ma as ma

archive1 = 'ich_ee.nc'
archive2 = 'ich_rk.nc'


nc1 = netCDF4.Dataset(archive1)
nc2 = netCDF4.Dataset(archive2)

time1 = nc1.variables['time'][:]

z1 = nc1.variables['depth'][:]
z2 = nc2.variables['depth'][:]

mortality1 = nc1.variables['mortality'][:]
mortality2 = nc2.variables['mortality'][:]

z1_m = z1.copy()
z2_m = z2.copy()

initime_ich=datetime.datetime(2006,1,1)
fintime_ich=[]
for t in time1:
    fintime_ich.append(initime_ich+datetime.timedelta(t/86400.))
    
for i in range(1, len(fintime_ich)):
    z1_m[i,:] = ma.masked_where(mortality1[i,:] == 4, z1_m[i,:])
    z2_m[i,:] = ma.masked_where(mortality2[i,:] == 4, z2_m[i,:])
    
medias1 = []
medias2 = []

for j in range(len(fintime_ich)):
    medias1.append(np.mean(z1_m[j,:]))
    medias2.append(np.mean(z2_m[j,:]))
    
medias1 = np.asarray(medias1)
medias2 = np.asarray(medias2)

sol = np.max(abs(medias1-medias2))
print(sol)
