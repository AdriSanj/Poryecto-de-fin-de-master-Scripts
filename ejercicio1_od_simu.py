# Basic imports for executing OpenDrift.
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from datetime import datetime, timedelta
import dateutil.parser

# For reading a .txt of initial positions if we want.
import pandas as pd

# Necessary for calculate the time 
import time


start = time.time()
filename = 'output.nc'
# Initial coordinates if we dont input a .txt file.
#longitude = -9
#latitude = 42

# Radius for uniform dispersion around source [meters]
R = 500.0


# Set number of floats 
N = 10000

# Time period of the simulation

start_time = datetime(2013, 4, 1)  
end_time = datetime(2013, 4, 16)



# Reading of the .txt file

names = ['lon', 'lat', 'dep']

df = pd.read_csv('particulas.txt', sep = ' ', names = names)

lon = df['lon']       
lat = df['lat']
z = -1*df['dep']



""" Model initialization """

model = OceanDrift(loglevel = 0)

""" Add readers """
# Ocean data

phys = reader_ROMS_native.Reader('modelo_hidrodinamico.nc')

# Adding readers
model.add_reader(phys)
model.set_config('general:use_auto_landmask', False)
object_type = 1 # PIW-1 Person-in-water (PIW), unknown state (mean values)

print('starting model run')
""" Seed elements """
model.seed_elements(lon = lon,
                lat = lat,
                number = N,
                #radius = R,
                #radius_type = "uniform",
                time = start_time,
                z = z      # comment if dont want depth
                #object_type=object_type
                )

# Para evitar el beaching
#model.set_config('general:coastline_action', 'previous')

model.set_config('drift:advection_scheme', 'runge-kutta4')
model.set_config('drift:horizontal_diffusivity', 53.14)


vm = False

if vm == True:
    model.set_config('drift:vertical_advection', True)
    model.set_config('drift:vertical_mixing', True)
    model.set_config('vertical_mixing:diffusivitymodel', 'environment')
    timestep_vm = 10     # segundos
    model.set_config('vertical_mixing:timestep', timestep_vm)
elif vm == False:
    model.set_config('drift:vertical_advection', False)
    model.set_config('drift:vertical_mixing', False)
# vertical mixing 
#model.set_config('drift:vertical_mixing', False)
#model.set_config('vertical_mixing:diffusivitymodel', 'environment')
#timestep_vm = 10     # segundos
#model.set_config('vertical_mixing:timestep', timestep_vm)

""" Model run """

# para pruebas, timestepdelta 30, modelo final 10 min
print('starting model run')
model.run(end_time = end_time, 
        time_step = timedelta(minutes = 10), 
        time_step_output = timedelta(minutes = 10), 
        outfile = filename,
        export_variables = ["trajectory", 
                            "time", 
                            "status",
                            "age_seconds",
                            "lon", 
                            "lat",
                            "z"]       
        )


print('finished model run')

# Exports the final time in which the program was executed.
end = time.time()
print('time of execution')
print(end-start)

f = open("tiempo_mask.txt", "w")
f.write(str(end-start))
f.close() 
