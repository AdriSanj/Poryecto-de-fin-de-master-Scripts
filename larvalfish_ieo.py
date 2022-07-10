# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2021, Trond Kristiansen, Niva
# Jan 2021 Simplified by Knut-Frode Dagestad, MET Norway, and adapted to to Kvile et al. (2018)

import datetime
import numpy as np
import logging; logger = logging.getLogger(__name__)
from opendrift.models.oceandrift_ieo import Lagrangian3DArray, OceanDrift


class LarvalFishElement(Lagrangian3DArray):
    """
    Extending Lagrangian3DArray with specific properties for larval and juvenile stages of fish
    """

    variables = Lagrangian3DArray.add_variables([
        ('diameter', {'dtype': np.float32,
                      'units': 'm',
                      'default': 0.0014}),  # for NEA Cod
        ('neutral_buoyancy_salinity', {'dtype': np.float32,
                                       'units': 'PSU',
                                       'default': 31.25}),  # for NEA Cod
        ('stage_fraction', {'dtype': np.float32,  # to track percentage of development time completed
                            'units': '',
                            'default': 0.}),
        ('hatched', {'dtype': np.uint8,  # 0 for eggs, 1 for larvae
                     'units': '',
                     'default': 0}),
        ('length', {'dtype': np.float32,
                    'units': 'mm',
                    'default': 2.79}),
        ('weight', {'dtype': np.float32,
                    'units': 'mg',
                    'default': 0.08}),
        ('egg_dens', {'dtype': np.float32,
                      'units': 'kgm-3',
                      'default': 1025}),
        ('dens_diff', {'dtype': np.float32,
                      'units': 'kgm-3',
                      'default': 0.}),
        ('zoopl', {'dtype': np.float32,
                      'units': 'mgm-3',
                      'default': 0.}),
        ('swim_speed', {'dtype': np.float32,
                      'units': '',
                      'default': 0.}),                  
        ('survival', {'dtype': np.float32,  # Not yet used
                      'units': '',
                      'default': 1.})])


class LarvalFish_sardine(OceanDrift):
    """Buoyant particle trajectory model based on the OpenDrift framework.

        Developed at MET Norway

        Generic module for particles that are subject to vertical turbulent
        mixing with the possibility for positive or negative buoyancy

        Particles could be e.g. oil droplets, plankton, or sediments

    """

    ElementType = LarvalFishElement

    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'sea_surface_wave_significant_height': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'zooplankton':{'fallback': 0},
        'land_binary_mask': {'fallback': None},
        'sea_floor_depth_below_sea_level': {'fallback': 100},
        'ocean_vertical_diffusivity': {'fallback': 0.01, 'profiles': True},
        'sea_water_temperature': {'fallback': 10, 'profiles': True},
        'sea_water_salinity': {'fallback': 34, 'profiles': True}
    }

    required_profiles_z_range = [0, -50]  # The depth range (in m) which profiles should cover

    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(LarvalFish_sardine, self).__init__(*args, **kwargs)

        # IBM configuration options
        self._add_config({
            'IBM:fraction_of_timestep_swimming':
                {'type': 'float', 'default': 0.15,
                 'min': 0.0, 'max': 1.0, 'units': 'fraction',
                 'description': 'Fraction of timestep swimming',
                 'level': self.CONFIG_LEVEL_ADVANCED},
            })

        self._set_config_default('drift:vertical_mixing', True)

    def update_terminal_velocity(self, Tprofiles=None,
                                 Sprofiles=None, z_index=None, DENSsp = 0.):

        g = 9.81  # ms-2
        eggs = np.where((self.elements.hatched == 0))[0] 
        if len(eggs) == 0:
            return

        # Pelagic Egg properties that determine buoyancy
        eggsize = self.elements.diameter  		
        eggsalinity = self.elements.neutral_buoyancy_salinity
        

        # prepare interpolation of temp, salt
        if not (Tprofiles is None and Sprofiles is None):
            if z_index is None:
                z_i = range(Tprofiles.shape[0])  # evtl. move out of loop
                # evtl. move out of loop
                z_index = interp1d(-self.environment_profiles['z'],
                                   z_i, bounds_error=False)
            zi = z_index(-self.elements.z)
            upper = np.maximum(np.floor(zi).astype(np.uint8), 0)
            lower = np.minimum(upper + 1, Tprofiles.shape[0] - 1)
            weight_upper = 1 - (zi - upper)

        # do interpolation of temp, salt if profiles were passed into
        # this function, if not, use reader by calling self.environment
        if Tprofiles is None:
            T0 = self.environment.sea_water_temperature 
            #if time_index == 0:
            #    T0_ini = self.environment.sea_water_temperature
        else:
            T0 = Tprofiles[upper, range(Tprofiles.shape[1])]*weight_upper + Tprofiles[lower, range(Tprofiles.shape[1])]*(1 - weight_upper)
            
        if Sprofiles is None:
            S0 = self.environment.sea_water_salinity
        else:
            S0 = Sprofiles[upper, range(Sprofiles.shape[1])]*weight_upper + Sprofiles[lower, range(Sprofiles.shape[1])]*(1 - weight_upper)

        # The density difference between a pelagic egg and the ambient water
        # is regulated by their salinity difference through the
        # equation of state for sea water.
        # The Egg has the same temperature as the ambient water and its
        # salinity is regulated by osmosis through the egg shell.
        
        # Modificacion. Calculo la densidad de spawning
	#ini_water_density = self.sea_water_density(T=T0ini,S=S0ini)
	#DENS_sp = 0.78797*ini_water_density+4.2519
        DENSw = self.sea_water_density(T=T0, S=S0)[eggs]
             
        self.elements.egg_dens[eggs] = self.dens_growth(self.elements.stage_fraction, DENSsp)[eggs]	# 
        print('STAGE FRACTION')
        print(self.elements.stage_fraction[0])
        #print(len(DENSegg))
        #print(DENSw)
        #print(DENSegg)
        dr = DENSw - self.elements.egg_dens[eggs]  # density difference
        self.elements.dens_diff[eggs] = dr
        
        # Velocidad acorde a Ichthyop
        
        #g_cm = 980        # Gravedad en cm/s2
        #mean_minor_axis = 0.05
        #mean_major_axis = 0.14
        #molec_visc = 0.01
        #LOGN = np.log(2 * mean_major_axis / mean_minor_axis)
        #W = (( g_cm * mean_minor_axis**2 / (24 * molec_visc * DENSw[eggs]) * (LOGN + 0.5) * dr )/100.)*self.time_step.total_seconds()

        # water viscosity
        my_w = 0.001 * (1.7915 - 0.0538 * T0[eggs] + 0.007 * (T0[eggs]** (2.0)) - 0.0023 * S0[eggs])
        # ~0.0014 kg m-1 s-1
        # terminal velocity for low Reynolds numbers
        W = (1.0 / my_w) * (1.0 / 18.0) * g * eggsize[eggs] ** 2 * dr

        # check if we are in a Reynolds regime where Re > 0.5
        highRe = np.where(W * 1000 * eggsize[eggs]/ my_w > 0.5)

        # Use empirical equations for terminal velocity in
        # high Reynolds numbers.
        # Empirical equations have length units in cm!
        my_w = 0.01854 * np.exp(-0.02783 * T0[eggs])  # in cm2/s
        d0 = (eggsize[eggs] * 100) - 0.4 * \
             (9.0 * my_w ** 2 / (100 * g) * DENSw / dr) ** (1.0 / 3.0)  # cm
        W2 = 19.0 * d0 * (0.001 * dr) ** (2.0 / 3.0) * (my_w * 0.001 * DENSw) ** (-1.0 / 3.0)
        # cm/s
        W2 = W2 / 100.  # back to m/s

        W[highRe] = W2[highRe]
        
        self.elements.terminal_velocity[eggs] = W


    def dens_growth(self, P_t, DENSsp):	
    	# Function to calculate the density of the egg with the time fraction and spawning density.
    	# Constants from eq. (4) Garcia et al.
        a_0 = 0.0012
        a_1 = -1.8528
        a_2 = 15.507
        a_3 = -41.1312 
        a_4 = 40.5289
        a_5 = -12.1166
        
        d_rho = 3.58 
        
        rho_n = a_0 + a_1*P_t + a_2*P_t**2 + a_3*P_t**3 + a_4*P_t**4 + a_5*P_t**5
        #print('dens_spawning')
        #print(DENSsp)
        #print(len(rho_n))
        
        rho_h = rho_n*d_rho + DENSsp	# Eq. (4) Garcia et al.
        #print('dens_growth')
        #print(rho_h[0])
        
        return rho_h


    def update_fish_larvae(self):
        A = 5389.61
        b = 1.59
        # Hatching of eggs
        eggs = np.where(self.elements.hatched==0)[0]
        #self.elements.length[eggs] = 2.8			# En mm
        #print(eggs)
        #print(len(eggs))
        if len(eggs) > 0:
            amb_duration = A*(self.environment.sea_water_temperature[eggs])**(-b) 	#CAMBIO
            #print(amb_duration)
            #print('temperatura del agua')
            #print(self.environment.sea_water_temperature[0])
            days_in_timestep = self.time_step.total_seconds()/(60.*60.)  # The fraction of a day completed in one time step.amd_duration viene dada en horas y asi anulamos unidades
            amb_fraction = days_in_timestep/amb_duration # Fraction of development time completed during present time step
            self.elements.stage_fraction[eggs] += amb_fraction # Add fraction completed during present timestep to cumulative fraction completed
            hatching = np.where(self.elements.stage_fraction[eggs]>=1)[0]
            if len(hatching) > 0:
                logger.debug('Hatching %s eggs' % len(hatching))
                self.elements.hatched[eggs[hatching]] = 1 # Eggs with total development time completed are hatched (1)

        larvae = np.where(self.elements.hatched>=1)[0]
        
        
        if len(larvae) == 0:
            logger.debug('%s eggs, with maximum stage_fraction of %s (1 gives hatching)'
                         % (len(eggs), self.elements.stage_fraction[eggs].max()))
            return
        
        larvae2 = np.where(self.elements.length[larvae] >= 4.5)[0]
        self.elements.hatched[np.where(self.elements.length >= 4.5)[0]] = 2
        if len(larvae) > 0 or len(larvae2) > 0:
            # Temperature in celsius acording to Lett et al. (2008)
            coef_1 = 1e-5
            coef_2 = 0.0308
            
            self.elements.length[np.where(self.elements.hatched==1)[0]] += (coef_1 + coef_2*self.environment.sea_water_temperature[np.where(self.elements.hatched == 1)[0]])*self.time_step.total_seconds()/86400.
            
            #Segundo estado de larva
            self.elements.zoopl = self.environment.zooplankton
            coeffz = 12.*8.			#12.*8.
            food = coeffz*self.elements.zoopl[np.where(self.elements.hatched == 2)[0]]
            self.elements.length[np.where(self.elements.hatched==2)[0]] += (coef_1 + coef_2*self.environment.sea_water_temperature[np.where(self.elements.hatched == 2)[0]])*food/(food + 0.1)*self.time_step.total_seconds()/86400.
         
         




    def update(self):

        self.update_fish_larvae()
        self.advect_ocean_current()
        if len(np.where(self.elements.hatched == 0)[0]) != 0:
            self.update_terminal_velocity()
        self.vertical_mixing()
        #self.larvae_vertical_migration()
