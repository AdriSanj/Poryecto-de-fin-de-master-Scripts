import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import datetime
import numpy.ma as ma

# ich

archive1 = 'sardine_ich.nc'
nc1 = netCDF4.Dataset(archive1)
z_ich = nc1.variables['depth'][:]
time1 = nc1.variables['time'][:]
mortality = nc1.variables['mortality'][:]
l1 = nc1.variables['length'][:]

#od

archive2 = 'sardine_od.nc'
nc2 = netCDF4.Dataset(archive2)
time2 = nc2.variables['time'][:]
z_od = nc2.variables['z'][:]
l2 = nc2.variables['length'][:]



initime_od=datetime.datetime(1970,1,1)
fintime_od=[]
for t in time2:
    fintime_od.append(initime_od+datetime.timedelta(t/86400.))

initime_ich=datetime.datetime(2006,1,1)
fintime_ich=[]
for t in time1:
    fintime_ich.append(initime_ich+datetime.timedelta(t/86400.))
    
z_od_m = z_od.copy()


for i in range(1, len(fintime_od)):
    z_od_m[:,i] = np.ma.filled(z_od_m[:,i].astype(float), z_od_m[:,i-1])   


z_ich_m = z_ich.copy()
l_ich_m = l1.copy()

for i in range(1, len(fintime_ich)):
    z_ich_m[i,:] = ma.masked_where(mortality[i,:] == 4, z_ich_m[i,:])
    l_ich_m[i,:] = ma.masked_where(mortality[i,:] == 4, l_ich_m[i,:])


# Comparacion particulas vivas

medias_OD = []
medias_ICH = []

mediasl_OD = []
mediasl_ICH = []


for i in range(len(fintime_od)):
    medias_OD.append(np.mean(z_od[:,i]))
    mediasl_OD.append(np.mean(l2[:,i]))
    
for j in range(len(fintime_ich)):
    medias_ICH.append(np.mean(z_ich_m[j,:]))
    mediasl_ICH.append(np.mean(l_ich_m[j,:]))


# Saltos de fase
hatch_1 = 306

hatch_2 = 755


# Profundidades medias

fig_1 = plt.figure()
plt.plot(fintime_od, medias_OD,'b-')
plt.plot(fintime_ich, medias_ICH,'r-')
plt.axvline(x = fintime_od[hatch_1], color = 'g')
plt.axvline(x = fintime_od[hatch_2], color = 'g')
plt.title('Profundidad vs Tiempo')
plt.xlabel('Tiempo')
plt.ylabel('Profundidad (m)')
plt.grid()
plt.show()    
fig_1.savefig('comp_depth_od_ich.png')
    
    
# Longitudes medias
    
fig_2 = plt.figure()
plt.plot(fintime_od, mediasl_OD,'b-')
plt.plot(fintime_ich, mediasl_ICH,'r-')
plt.axvline(x = fintime_od[hatch_1], color = 'g')
plt.axvline(x = fintime_od[hatch_2], color = 'g')
plt.title('Longitud vs Tiempo')
plt.xlabel('Tiempo')
plt.ylabel('Longitud (mm)')
plt.grid()
plt.show()

fig_2.savefig('comp_length_od_ich.png')
