##################################################################
#                                                                #
# tSPICE: Tidal Signal with Python and SPICE                     #
#                                                                #
##################################################################
# License: GNU Affero General Public License v3 (AGPL-3.0)       #
# Authors: Jorge I. Zuluaga, Juanita A. Agudelo                  #
# Contact: jorge.zuluaga@udea.edu.co                             #
##################################################################
"""
Tidal Signal Calculations Module.

This module provides classes and methods for calculating tidal potentials and
internal deformations of celestial bodies using SPICE ephemerides.

Classes
-------
Body
    Represents a celestial body for tidal potential calculations.
BodyResponse
    Extends Body with internal structure modeling and Love number calculations.

Examples
--------
Calculate tidal potential from the Moon and Sun on Earth:

>>> import tspice
>>> tspice.initialize()
>>> earth = tspice.Body('Earth')
>>> loc = {'lat': 4.49, 'lon': -73.14, 'depth': 0}
>>> dates = {'start': '2025-01-01', 'stop': '2025-01-02', 'step': '1h'}
>>> tgp = earth.tgp_many_bodies(['Moon', 'Sun'], loc, dates)

Calculate Love numbers for Earth:

>>> from tspice import BodyResponse, Earth
>>> earth_response = BodyResponse('Earth')
>>> earth_model = Earth()
>>> earth_response.set_integration_parameters_ad(
...     n=2, f_days=1.0,
...     layers_list=[...], 
...     planet_profile=earth_model.planet_profile
... )
>>> earth_response.integrate_internal_solutions_ad()
>>> print(f"Love number h_2 = {earth_response.h_n}")
"""

# Import modules
from tspice.utils import *

# Import libraries
import numpy as np
import spiceypy as spy
import scipy.special as sps
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os

# Compatibility for new SciPy versions (e.g. >= 1.15.0) where sph_harm is removed
if not hasattr(sps, 'sph_harm'):
    if hasattr(sps, 'sph_harm_y'):
         # Map old signature sph_harm(m, n, theta=azimuthal, phi=polar)
         # to new signature sph_harm_y(n, m, theta=polar, phi=azimuthal)
         def _sph_harm_compat(m, n, theta, phi):
             return sps.sph_harm_y(n, m, phi, theta)
         sps.sph_harm = _sph_harm_compat

# File paths
file_dir = os.path.dirname(os.path.abspath(__file__)) #Current file location
data_dir = os.path.join(file_dir, 'data')	#Data files

class Body():
	def __init__(self, name, solar_system=True):
		
		'''
		This class defines the Body for tidal calculations. The body could be a Solar System body or an exoplanet from a real or hypothetical system.
            
		Input:
		- name: Name of the body, e.g. 'Earth', 'Mars', 'Io', 'Wasp-12b', etc.
		- solar_system: Boolean to indicate if the body is in the solar system and made part of the SPICE kernels (default True) or an exoplanet from a real or hypothetical system (False). In the case of exoplanets, the user should provide the necessary physical and orbital parameters later for the tidal calculations.
        
        Output:
		- Body object with the physical parameters of the body.
		'''

		#By default, the bodies are from the solar system
		if solar_system:
		
		#Check if kernels are loaded before allowing creation of a Body
			import tspice
			if not tspice.kernels_loaded:
				raise RuntimeError(
					"SPICE kernels not loaded. Please call tspice.initialize() first. "
					"Example:\n\n"
					"import tspice\n"
					"tspice.initialize()\n"
					"# or specify a custom directory: tspice.initialize('/path/to/my/kernels')\n"
					"body = tspice.Body('Moon')"
				)
			
			#Name, Id and body-fixed frame of the (target) body
			self.body_name = name
			self.body_id = spy.bodn2c(name)
			self.body_frame = spy.cidfrm(self.body_id)[1]
			
			try:
				#GM for the body
				self.GM = spy.bodvrd(name, 'GM', 1)[1][0]    #in [km^3/s^2]
			except Exception as e:
				print(f"Error getting GM for {name}: {e}")
				self.GM = None
				
			try:
				#a_p is the mean radius for the body
				self.a_ellips = spy.bodvrd(name, 'RADII', 3)[1]		#in [km]
				self.a_p = self.a_ellips.mean()	#in [km]
				#Flattening of the body
				self.f = (self.a_ellips[0]-self.a_ellips[-1])/self.a_ellips[0]
			except Exception as e:
				print(f"Error getting RADII for {name}: {e}")
				self.a_p = None

			if self.GM and self.a_p:
				self.g_ref = 1e3*self.GM/self.a_p**2	#in [m/s^2]
		
		else:
			#[PENDING] Add functionality for other bodies outside the solar system
			raise NotImplementedError("Currently, only solar system bodies are supported.")

			#[PENDING] Add physical parameters for exoplanets or hypothetical bodies

	def array_et_utc(self, date, return_step_seconds=False):

		'''
		This function creates an array of times in ET from a dictionary with start, stop, step, and time_frame (UTC or TDB).

		Input:
		- date: Dictionary with start, stop, step, and time_frame (UTC or TDB).

		Output:
		- et_utc: Array of times in ET.
		'''

		if date['time_frame'] == 'UTC' or 'time_frame' not in date.keys():
			et_utc_start = spy.utc2et(date['start'])
			et_utc_stop = spy.utc2et(date['stop'])
			step_seconds = convert_step_to_seconds(date['step'])
			et_utc = np.arange(et_utc_start, et_utc_stop, step_seconds)

		#[PENDING] Add more conditions later
		else:
			raise ValueError("Currently, only 'UTC' time frame is supported.")

		if return_step_seconds:
			return step_seconds, et_utc
		else:
			return et_utc

	#[PENDING] Function to get the body IDs from their names

	def Kn_func(self, n, M_ext, M=None, a_p=None):

		'''
		This program calculate the constant Kn in the Tide-Generating Potential expansion.
		
		Input:
		- n: Degree of the term
		- a_p: Radius of the target body [m]
		- M: Mass of the target body [any mass units]
		- M_ext: Mass of the external body [the same as M]

		Output:
		- Kn: Constant K_n for a = a_p
		
		Note:
		1. To get the potential for other distance from the COM,
		should multiply this value by (a/a_p)^n
		2. This definition is different to the one use in Agnew (2015)
		that includes the factor 1/r_mean, with r_mean the mean distance
		to the external body.'''

		if M is None and a_p is None:
			M = self.GM
			a_p = self.a_p

		Kn = (4*np.pi/(2*n+1))*(M_ext/M)*a_p**(n+2)	#in [km^(n+2)]
		return Kn
	
	def subpoint_coordinates(self, et_utc, body):

		'''
		This function calculates the geographical coordinates of the subpoint
		of an external body on the target body for an array of times in ET.
		
		Input:
		- et_utc: Array of times in ET.
		- body: Name of the external body, e.g. 'Moon', 'Sun', 'Mercury', etc.

		Output:
		- phis_ext: Array of longitudes of the subpoint [rad].
		- thetas_ext: Array of colatitudes of the subpoint [rad].
		- alts_ext: Array of altitudes of the subpoint above the reference ellipsoid [km].
		'''

		#Loop over the times to get the coordinates of the external body subpoint
		lons_ext, lats_ext, alts_ext = np.zeros(len(et_utc)), np.zeros(len(et_utc)), np.zeros(len(et_utc))
		for i, t in enumerate(et_utc):

			#Rectangular coordinates of the subpoint
			subp, et, pos = spy.subpnt('INTERCEPT/ELLIPSOID', self.body_name, t, self.body_frame, 'XCN+S', body)

			#Geographical coordinates of the subpoint
			lon, lat, alt = spy.recgeo(subp, self.a_p, self.f) #in [rad, rad, km]
			lons_ext[i], lats_ext[i], alts_ext[i] = lon, lat, alt
		
		#Final coordinates arrays
		phis_ext, thetas_ext = lons_ext, np.pi/2 - lats_ext

		return phis_ext, thetas_ext, alts_ext

	def tgp_one_body(self, body, loc_sta=None, dates=None, nmax=None, time_array=False):

		'''
		This function calculates the TGP for a point on the target body with geographic coordinates (lon_s, lat_s) at a distance a from the COM,
		due to a list of external bodies.
		
		Input:
		- loc_sta: Dictionary with the geographic coordinates of the station, including depth. For example: dict(lat=4.49, lon=-73.14, depth=0)
		- dates: Dictionary with start, stop, and step dates. For example: dict(start="2025-01-01 00:00:00", stop="2025-01-29 00:00:00", step="1h")
		- bodies: the body to consider, e.g. 'Moon', 'Sun', 'Mercury'.
		- nmax: Maximum degree to consider in the tidal calculation.

		Returns:
		- Vtid: Array with the time series of V/g at the station for each body
		- tjd: Times in JD of the ephemerides (and of V/g).
		'''

		#Cordinates of the station
		if loc_sta is None:
			#If the function is called within tgp_many_bodies
			phi_sta, theta_sta, a_sta = self.phi_sta, self.theta_sta, self.a_sta
		else:
			phi_sta, theta_sta, a_sta = loc_func(loc_sta, self.a_ellips)

		#Array of UTC times in ET
		if dates is None:
			#If the function is called within tgp_many_bodies
			et_utc = self.et_utc
		else:
			et_utc = self.array_et_utc(dates)

		#Maximum degree
		if nmax is None:
			nmax = self.nmax
		else:
			nmax = nmax

		#GM for the external body
		GM_ext = spy.bodvrd(body, 'GM', 1)[1][0] #in [km^3/s^2]

		#State vector and light time travel
		x, lt = spy.spkezr(body, et_utc, 'ECLIPJ2000', 'XCN+S', self.body_name)	#in [km]

		#Position vector
		r_ext = np.array(x)[:,:3]	#in [km]
		r = np.linalg.norm(r_ext, axis=1)
		#r_mean = r.mean()
		#xit = r/r_mean	#Distance variation in time

		#Subpoint coordinates of the external body
		phis_ext, thetas_ext, alts_ext = self.subpoint_coordinates(et_utc, body)

		#Vtid/g for all the times
		V_g_tid_body = np.zeros(len(et_utc))

		#Sum over the degrees
		n_init = 2
		n_final = nmax
		for n in range(n_init,n_final+1):
			Kn = self.Kn_func(n, GM_ext, self.GM, self.a_p)
			aamain_n = (a_sta/self.a_p)**n	#For a != a_p
			fn = aamain_n*Kn/r**(n+1)	#in [km] because a and r is in [km]

			#Sum over the orders
			sigma_n = 0
			for m in range(-n,n+1):
					Ynm_ext = sps.sph_harm(m, n, phis_ext, thetas_ext)  #External body
					Ynm_sta = sps.sph_harm(m, n, phi_sta, theta_sta)  #Earth (station)
					YeYs = np.conjugate(Ynm_ext)*Ynm_sta
					sigma_n += YeYs.real

			#Potential contribution of degree n
			V_g_tid_body += fn*sigma_n	
		
		#Factor 1e3 to return in [m]
		V_g_tid_body = V_g_tid_body*1e3

		if time_array:
			return V_g_tid_body, et_utc
		else:
			return V_g_tid_body
	
	def tgp_many_bodies(self, bodies, loc_sta, dates, nmax=6, time_array=False, body_signal=False, verbose=True):

		'''
		This function calculates the TGP for a point on the target body with geographic coordinates (lon_s, lat_s) at a distance a from the COM, due to a list of external bodies.
		
		Input:
		- loc_sta: Dictionary with the geographic coordinates of the station, including depth. For example: dict(lat=4.49, lon=-73.14, depth=0)
		- dates: Dictionary with start, stop, and step dates. For example: dict(start="2025-01-01 00:00:00", stop="2025-01-29 00:00:00", step="1h")
		- bodies: List of external bodies to consider, e.g. ['Moon', 'Sun', 'Mercury']
		- nmax: Maximum degree to consider in the tidal calculation.

		Returns:
		- Vtid: Array with the time series of V/g at the station for each body
		- tjd: Times in JD of the ephemerides (and of V/g).
		'''

		#Planetary ID for external bodies (outer planets here are really the barycenters)
		#bodies_ids = {'Moon':301, 'Sun':10, 'Mercury':199,
					#'Venus':299, 'Earth':399, 'Mars':4, 'Jupiter':5,
					#'Saturn':6, 'Uranus':7, 'Neptune':8}

		#Cordinates of the station
		phi_sta, theta_sta, a_sta = loc_func(loc_sta, self.a_ellips)
		self.phi_sta, self.theta_sta, self.a_sta = phi_sta, theta_sta, a_sta

		#Array of UTC times in ET
		et_utc = self.array_et_utc(dates)
		self.et_utc = et_utc

		#Maximum degree
		self.nmax = nmax

		#Array to store all V/g time series for all bodies
		V_g_tid_array = np.zeros((len(et_utc),len(bodies)))

		#Loop over external bodies
		for i,body in enumerate(bodies):
			#V/g time series for the body
			V_g_tid_body = self.tgp_one_body(body,)
			V_g_tid_array[:,i] = V_g_tid_body
			if verbose: print(f'{body} contribution calculated!')
		
		if body_signal == True:
			if time_array:
				return V_g_tid_array, et_utc
			else:
				return V_g_tid_array
		else:
			#Total V/g summing the contributions of all bodies
			V_g_tid_total = V_g_tid_array.sum(axis=1)
			if time_array:
				return V_g_tid_total, et_utc
			else:
				return V_g_tid_total

	def plot_one_signal(self, et, tgp, loc, colors, label=r'V_{\text{tid}}(t)/g', units='cm', mean_value=False, savepath=None):
		'''
		Plot a single tidal signal over time.
		
		Input:
		- et: Array of ephemeris times.
		- tgp: Tidal signal array.
		- loc: Dictionary with location (lat, lon, depth).
		- colors: List of colors for the plot.
		- label: Label for y-axis (default: r'V_{\text{tid}}(t)/g').
		- units: Units for the signal (default: 'cm').
		- mean_value: Whether to plot mean value line (default: False).
		- savepath: Path to save the figure (default: None).
		'''
		utc_times = et_to_utc_string(et)
		ticks = np.linspace(0, len(utc_times)-1, 6, dtype=int)

		fig, ax = plt.subplots(figsize=(8,4))
		ax.plot(utc_times, tgp, color=colors[0])

		ax.set_xlabel('Time [UTC]', fontsize=12)
		ax.set_ylabel(r'$'+label+'$ ['+units+']', fontsize=12)
		ax.set_title('Lon $= %.2f$°; Lat $= %.2f$°; Depth $= %.2f$ km'%(loc['lon'], loc['lat'], loc['depth']), fontsize=12)
		ax.set_xticks([utc_times[i] for i in ticks])
		ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
		ax.xaxis.set_tick_params(direction='in', which='both', labelsize=12)
		ax.yaxis.set_tick_params(direction='in', which='both', labelsize=12)
		ax.grid(alpha=0.2)
		ax.margins(x=0)

		if mean_value:
			tgp_mean = tgp.mean()
			ax.hlines(tgp_mean, utc_times[0], utc_times[-1], color=colors[1], label=r'Mean$=%.2f$ %s'%(tgp_mean, units))
			ax.legend(fontsize=12, loc='upper right')

		fig.tight_layout()

		if savepath:
			fig.savefig(savepath, bbox_inches='tight', dpi=300)

	def plot_many_signal(self, et, tgps, loc, colors, y_label=r'V(t)/g', signal_labels=None, units='cm', savepath=None):
		'''
		Plot multiple tidal signals on the same plot.
		
		Input:
		- et: Array of ephemeris times.
		- tgps: List of tidal signal arrays.
		- loc: Dictionary with location (lat, lon, depth).
		- colors: List of colors for each signal.
		- y_label: Label for y-axis (default: r'V(t)/g').
		- signal_labels: List of labels for each signal (default: None).
		- units: Units for the signals (default: 'cm').
		- savepath: Path to save the figure (default: None).
		'''
		utc_times = et_to_utc_string(et)
		ticks = np.linspace(0, len(utc_times)-1, 6, dtype=int)

		fig, ax = plt.subplots(figsize=(8,4))
		for i,signal in enumerate(tgps):
			ax.plot(utc_times, signal, color=colors[i], label=signal_labels[i] if signal_labels else None)

		ax.set_xlabel('Time [UTC]', fontsize=12)
		ax.set_ylabel(r'$'+y_label+'$ ['+units+']', fontsize=12)
		ax.set_title('Lon $= %.2f$°; Lat $= %.2f$°; Depth $= %.2f$ km'%(loc['lon'], loc['lat'], loc['depth']), fontsize=12)
		ax.set_xticks([utc_times[i] for i in ticks])
		ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
		ax.xaxis.set_tick_params(direction='in', which='both', labelsize=12)
		ax.yaxis.set_tick_params(direction='in', which='both', labelsize=12)
		ax.grid(alpha=0.2)
		ax.margins(x=0)
		if signal_labels:
			ax.legend(fontsize=12, loc='best')

		fig.tight_layout()

		if savepath:
			fig.savefig(savepath, bbox_inches='tight', dpi=300)

class BodyResponse(Body):

    def __init__(self, name):

        #The class inherit the attributes from Body
        super().__init__(name)
        #GM_main, a_ellips, a_main, f_main, g_ref

    def scale_constants(self, verbose=True):

        '''
        Function to calculate characteristic scales for adimensionalization of the internal solutions yis.

        Inputs:
        - verbose: Boolean to print the characteristic scales (default: True)

        Outputs:
        - L: Length scale [m]
        - M: Mass scale [kg]
        - RHO: Density scale [kg/m^3]
        - P: Pressure/Elasticity scale [Pa]
        - V: Velocity scale [m/s]
        - T: Time scale [s]
        - OMEGA: Angular frequency scale [rad/s]
        - Gad: Gravity scale [m/s^2]

        '''

        #Scale constants
        self.G = 6.67430e-11	#Gravitational constant in [m^3/kg/s^2]

        #[PENDING] Check this, I'm making int(self.a_main) for the radius to be the same as the PREM model
        #In the future we should guarantee that the planetary model has the same radius as the SPICE mean radius
        self.L = int(self.a_p)*1e3	#Length in [m] 

        self.M = self.GM*1e9/self.G	#Mass in [kg]
        self.RHO = self.M/self.L**3	#kg/m^3
        self.P = self.G*self.M**2/self.L**4	#Pressure/Elasticity modules in [Pa]
        self.V = (self.P/self.RHO)**0.5	#Velocity in [m/s]
        self.T = self.L/self.V	#Time in [s]
        self.OMEGA = 1/self.T	#Angular frecuency in [rad/s]
        self.Gad = self.G*self.RHO*self.L #Gravity acceleration in [m/s^2]

        if verbose:
            print('Characteristic scales for adimensionalization:')
            print(f'Length scale L = {self.L:.2e} m')
            print(f'Mass scale M = {self.M:.2e} kg')
            print(f'Density scale RHO = {self.RHO:.2e} kg/m^3')
            print(f'Pressure/Elasticity scale P = {self.P:.2e} Pa')
            print(f'Velocity scale V = {self.V:.2e} m/s')
            print(f'Time scale T = {self.T:.2e} s')
            print(f'Angular frequency scale OMEGA = {self.OMEGA:.2e} rad/s')
            print(f'Gravity scale Gad = {self.Gad:.2e} m/s^2')


    def read_layers(self, layers_list, r0_ini_ad, nsteps_total=5e4, dr=1):

        '''
        Function to read the layers of the planetary model and define the steps for each layer in the integration.

        Inputs:
        - layers_list: List of layers for the integration. Each layer is a dictionary with the following keys:
            - 'name': Name of the layer (string)
            - 'r0': Initial radius of the layer [m]
            - 'rf': Final radius of the layer [m]
            - 'type': Type of the layer ('solid' or 'liquid')
        - r0_ini_ad: Initial radius for integration [adimensional]. Default is 1 km to avoid singularities at the center
        - nsteps_total: Number of steps for the total integration (default: 5e4 steps)
        - dr: Delta r to avoid discontinuities between layers [m]. Default is 1 m

        Outputs:
        - sorted_layers_list: List of layers sorted by initial radius (from center to surface)
        - nsteps_dict: Dictionary with the number of steps for each layer
        - steps_dict: Dictionary with the integration steps for each layer
        '''

        #Sorted layers by initial radius
        sorted_layers_list = sorted(layers_list, key=lambda x: x['r0'])

        #For the transition between layers we choose a step of 1 meter in adimensional units
        #dr_ad = dr/self.L

        #Calculate the steps for each layer in the integration
        nsteps_dict = {}
        steps_dict = {}
        for l,layer in enumerate(sorted_layers_list):

            #For the first layer
            if l==0:
                rf_ad = layer['rf']/self.L	#[adimensional]
                nsteps = int(rf_ad*nsteps_total)
                nsteps_dict[layer['name']] = nsteps
                steps_dict[layer['name']] = np.linspace(r0_ini_ad, rf_ad, nsteps)

            #For the other layers
            else:
                r0_ad = (layer['r0']+dr)/self.L	#[adimensional]
                rf_ad = layer['rf']/self.L	#[adimensional]
                nsteps = int((rf_ad - r0_ad)*nsteps_total)
                nsteps_dict[layer['name']] = nsteps
                steps_dict[layer['name']] = np.linspace(r0_ad, rf_ad, nsteps)

        return sorted_layers_list, nsteps_dict, steps_dict
    
    def set_integration_parameters_ad(self, n, f_days, layers_list, planet_profile, nsteps=5e4, r0_ini=1e3):

        '''
        Define the adimensional parameters needed for the integration of the equations of motion: n, omega, planetary profiles, steps, layers, etc.
        
        Inputs:
        - n: Degree of the tidal solution
        - f_days: Frequency of the tidal forcing [cycles per day]
        - nsteps: Number of steps for the total integration
        - r0_ini: Initial radius for integration [m]. Default is 1 km to avoid singularities at the center
        - layers_list: List of layers for the integration. Each layer is a dictionary with the following keys:
            - 'name': Name of the layer (string)
            - 'r0': Initial radius of the layer [m]
            - 'rf': Final radius of the layer [m]
            - 'type': Type of the layer ('solid' or 'liquid')
        - planet_profile: Dictionary with the planetary profiles as functions of the radius (can be dimensional or adimensional). The dictionary must contain the following keys:
            - 'dimensionless': Boolean indicating if the profiles are adimensional (True) or dimensional (False)
            - 'rho': Density profile function [kg/m^3 or adimensional]
            - 'g': Gravitational acceleration profile function [m/s^2 or adimensional]
            - 'mu': Shear modulus profile function [Pa or adimensional]
            - 'lamb': Lamé's first parameter profile function [Pa or adimensional]

        Outputs:
        - None (the parameters are stored as attributes of the class). We can access them as self.n, self.omega, self.rho0_ad, etc.

        '''

        #Calculate the scale constants if not done yet  
        self.scale_constants(verbose=False)      
        
        #Degree of the solution
        self.n = n

        #Frequency of the tidal forcing
        self.f_days = f_days	#[cycles per day]
        self.omega = f_days*2*np.pi/(3600*24)	#[rad/s]
        self.omega_ad = self.omega/self.OMEGA	#Adimensional angular frequency

        #Planetary profiles (functions of r_ad)
        if planet_profile['dimensionless']:    #If the profiles are already adimensional
            self.rho0_ad = planet_profile['rho']
            self.mu_ad = planet_profile['mu']
            self.lamb_ad = planet_profile['lamb']
            self.g0_ad = planet_profile['g']
        else:   #If the profiles are in dimensional units we need to reescale them
            self.rho0_ad = lambda r_ad: planet_profile['rho'](r_ad*self.L)/self.RHO
            self.mu_ad = lambda r_ad: planet_profile['mu'](r_ad*self.L)/self.P
            self.lamb_ad = lambda r_ad: planet_profile['lamb'](r_ad*self.L)/self.P
            self.g0_ad = lambda r_ad: planet_profile['g'](r_ad*self.L)/self.Gad

        #Number of layers in our planetary model
        self.nsteps_total = nsteps

        #Initial radius for integration
        self.r0_ini = r0_ini    #[m]
        self.r0_ini_ad = self.r0_ini/self.L	#[adimensional]

        #Steps for each layer in the integration
        self.layers_list, self.nsteps_dict, self.steps_dict = self.read_layers(layers_list, self.r0_ini_ad, self.nsteps_total)

        #Adimensional parameters for the integration
        self.params_ad = {
            'n': self.n,
            'omega': self.omega_ad,
            'rho': self.rho0_ad,
            'lam': self.lamb_ad,
            'mu': self.mu_ad,
            'g': self.g0_ad
        }

    #[PENDING] Update integration parameters

    def initial_conditions_ad(self, r0_ini_ad, params_ad, setup='AmorinGudkova2024'):
        
        '''
        Function to set the initial conditions of the integration at r = r0, for three independent solutions.
        
        Inputs:
        - r0_ini_ad: Initial radius from the center to start the integration [adimensional]
        '''

        #Get the initial conditions based on the selected setup
        if setup == 'TakeuchiSaito1972':
            Y0_ad = Y0_TakeuchiSaito1972_ad(r0_ini_ad, params_ad)
        elif setup == 'XuSun2003':
            Y0_ad = Y0_XuSun2003_ad(r0_ini_ad, params_ad)
        elif setup == 'AmorinGudkova2024':
            Y0_ad = Y0_AmorinGudkova2024_ad(r0_ini_ad, params_ad)
        return Y0_ad

    def integrate_internal_solutions_ad(self, setup='AmorinGudkova2024', verbose=False):

        '''
        Function to integrate the equations of motion for a planet
        '''

        #Degree of the solution
        n = self.n

        #Initial conditions at r0
        self.Y0 = self.initial_conditions_ad(self.r0_ini_ad, self.params_ad, setup=setup)
        Y0_1_ad, Y0_2_ad, Y0_3_ad = self.Y0[:,0], self.Y0[:,1], self.Y0[:,2]

        #Integrate the equations of motion
        if setup == 'AmorinGudkova2024':

            if len(self.layers_list) == 3:
                pass  #Continue with the integration
            else:
                raise ValueError('For the Amorin & Gudkova (2024) setup we only support three layers (solid-fluid-solid). An additional (solid) layer can be added to the last one.')
            
            ###Integration through the layers###

            #First solid layer (three independent solutions)
            if self.layers_list[0]['type'] == 'solid':

                name = self.layers_list[0]['name']
                rs_ad = self.steps_dict[name]

                #Integration in the solid part
                y_1_inner_ad = solve_ivp(dydr_solid_AmorinGudkova2024_ad, (rs_ad[0], rs_ad[-1]), Y0_1_ad, t_eval=rs_ad, method='BDF', args=(self.params_ad,), dense_output=False)
                y_2_inner_ad = solve_ivp(dydr_solid_AmorinGudkova2024_ad, (rs_ad[0], rs_ad[-1]), Y0_2_ad, t_eval=rs_ad, method='BDF', args=(self.params_ad,), dense_output=False)
                y_3_inner_ad = solve_ivp(dydr_solid_AmorinGudkova2024_ad, (rs_ad[0], rs_ad[-1]), Y0_3_ad, t_eval=rs_ad, method='BDF', args=(self.params_ad,), dense_output=False)

                if verbose: print(f'Integrated first solid layer!')
            
            else:
                raise ValueError('The first layer must be solid for the Amorin & Gudkova (2024) setup.')
            
            #Second fluid layer (combine the three solutions)
            if self.layers_list[1]['type'] == 'fluid':

                name = self.layers_list[1]['name']
                rs_ad = self.steps_dict[name]

                #Gravity, density and radius at the interface solid-liquid (ICB)
                rs_icb = self.steps_dict[self.layers_list[0]['name']][-1]
                g0_icb = self.params_ad['g'](rs_icb)
                rho0_icb = self.params_ad['rho'](rs_icb)

                #We use the previous yis to define the new zis
                y1_A_c, y2_A_c, y4_A_c, y5_A_c, y6_A_c = y_1_inner_ad.y[[0,1,3,4,5],-1]
                y1_B_c, y2_B_c, y4_B_c, y5_B_c, y6_B_c = y_2_inner_ad.y[[0,1,3,4,5],-1]
                y1_C_c, y2_C_c, y4_C_c, y5_C_c, y6_C_c = y_3_inner_ad.y[[0,1,3,4,5],-1]

                #New coeficients to combine the three previous solutions
                denom = g0_icb*rho0_icb*(y1_B_c*y4_C_c - y4_B_c*y1_C_c) + (y4_B_c*y2_C_c - y2_B_c*y4_C_c) + rho0_icb*(y4_B_c*y5_C_c - y5_B_c*y4_C_c)
                numV = g0_icb*rho0_icb*(y4_A_c*y1_C_c - y1_A_c*y4_C_c) + (y2_A_c*y4_C_c - y4_A_c*y2_C_c) + rho0_icb*(y5_A_c*y4_C_c - y4_A_c*y5_C_c)
                numS = g0_icb*rho0_icb*(y1_A_c*y4_B_c - y4_A_c*y1_B_c) + (y4_A_c*y2_B_c - y2_A_c*y4_B_c) + rho0_icb*(y4_A_c*y5_B_c - y5_A_c*y4_B_c)
                V_A = numV/denom
                S_A = numS/denom

                #Functions to combine the yis then
                B_from_A = lambda A: A*V_A
                C_from_A = lambda A: A*S_A

                #New variables (combinations of the yis)
                z1_til = y1_A_c + V_A*y1_B_c + S_A*y1_C_c
                z5_til = y5_A_c + V_A*y5_B_c + S_A*y5_C_c
                z6_til = y6_A_c + V_A*y6_B_c + S_A*y6_C_c
                z7_til = z6_til + 4*np.pi*rho0_icb*z1_til + ((n+1)/rs_icb - 4*np.pi*rho0_icb/g0_icb)*z5_til

                #In fluid part, we just do one integration
                Z0_outer_ad = np.array([z5_til, z7_til])	

                #Integration in the fluid layer
                z_outer_ad_redu = solve_ivp(dzdr_fluid_AmorinGudkova2024_ad, (rs_ad[0], rs_ad[-1]), Z0_outer_ad, t_eval=rs_ad, method='BDF', args=(self.params_ad,), dense_output=False)

                if verbose: print(f'Integrated second fluid layer!')
            
            else:
                raise ValueError('The second layer must be fluid for the Amorin & Gudkova (2024) setup.')

            #Third solid layer (we have three unknown coefficients and solutions)
            if self.layers_list[2]['type'] == 'solid':

                name = self.layers_list[2]['name']
                rs_ad = self.steps_dict[name]

                #Gravity, density and radius at the interface liquid-solid (CMB)
                rs_cmb = self.steps_dict[self.layers_list[1]['name']][-1]
                g0_cmb = self.params_ad['g'](rs_cmb)
                rho0_cmb = self.params_ad['rho'](rs_cmb)

                #We use the previous zis to define the new yis
                z5_b, z7_b = z_outer_ad_redu.y[0,-1], z_outer_ad_redu.y[1,-1]
                
                #New three independent solutions from the zis
                Y0_alpha_mantle_ad = np.array([0, -rho0_cmb*z5_b, 0, 0, z5_b, z7_b - ((n+1)/rs_cmb - 4*np.pi*rho0_cmb/g0_cmb)*z5_b])	#Coefficient A
                Y0_beta_mantle_ad = np.array([1, rho0_cmb*g0_cmb, 0, 0, 0, -4*np.pi*rho0_cmb])	#Coefficient D=z1_b
                Y0_gamma_mantle_ad = np.array([0, 0, 1, 0, 0, 0])	#Coefficient E=y3_b / discontinuity in the tangential displacement

                #Integration in the solid part
                y_alpha_mantle_ad = solve_ivp(dydr_solid_AmorinGudkova2024_ad, (rs_ad[0], rs_ad[-1]), Y0_alpha_mantle_ad, t_eval=rs_ad, method='BDF', args=(self.params_ad,), dense_output=True)
                y_beta_mantle_ad = solve_ivp(dydr_solid_AmorinGudkova2024_ad, (rs_ad[0], rs_ad[-1]), Y0_beta_mantle_ad, t_eval=rs_ad, method='BDF', args=(self.params_ad,), dense_output=True)
                y_gamma_mantle_ad = solve_ivp(dydr_solid_AmorinGudkova2024_ad, (rs_ad[0], rs_ad[-1]), Y0_gamma_mantle_ad, t_eval=rs_ad, method='BDF', args=(self.params_ad,), dense_output=True)

                if verbose: print(f'Integrated third solid layer!')

            else:
                raise ValueError('The third layer must be solid for the Amorin & Gudkova (2024) setup.')
            
            ###Combining solutions###
            n_cond = 5  #Number of conditions

            #Matrix to store the conditions and system of equations
            Pmat = np.zeros((n_cond,n_cond))
            Bmat = np.zeros(n_cond)

            #Conditions:

            #y1(a) definition
            Pmat[0,:] = y_alpha_mantle_ad.y[0,-1], y_beta_mantle_ad.y[0,-1], y_gamma_mantle_ad.y[0,-1], -1, 0

            #Null stresses at surface
            Pmat[1,:] = y_alpha_mantle_ad.y[1,-1], y_beta_mantle_ad.y[1,-1], y_gamma_mantle_ad.y[1,-1], 0, 0 #y2(a)=0
            Pmat[2,:] = y_alpha_mantle_ad.y[3,-1], y_beta_mantle_ad.y[3,-1], y_gamma_mantle_ad.y[3,-1], 0, 0 # y4(a)=0

            #y5(a) definition
            Pmat[3,:] = y_alpha_mantle_ad.y[4,-1], y_beta_mantle_ad.y[4,-1], y_gamma_mantle_ad.y[4,-1], 0, -1

            #Potential continuity at surface
            Pmat[4,:] = y_alpha_mantle_ad.y[5,-1], y_beta_mantle_ad.y[5,-1], y_gamma_mantle_ad.y[5,-1], 0, (n+1)
            Bmat[4] = (2*n+1)

            #Solving the system of equations
            Cmat = np.linalg.solve(Pmat, Bmat)  #Coefficients
            A, D, E, y1_a_ad, y5_a_ad = Cmat 

            #Concatenate all the results
            rs_all_ad = [self.steps_dict[layer['name']] for layer in self.layers_list]
            self.rs_all_ad = np.concatenate(rs_all_ad)

            #Combine solutions:
            
            #In the inner core
            y_inner_solution_ad = A*y_1_inner_ad.y + B_from_A(A)*y_2_inner_ad.y + C_from_A(A)*y_3_inner_ad.y

            #In the outer core
            y_outer_solution_ad = np.zeros((6, self.nsteps_dict[self.layers_list[1]['name']]))   #y1, y2, y3, y6 are not computed in the outer core for this setup
            y_outer_solution_ad[3,:] = 0    #y4=0 in fluid outer core
            y_outer_solution_ad[4,:] = A*z_outer_ad_redu.y[0,:]     #y5=A*z5_outer_ad

            #Combine solutions in the mantle
            y_mantle_solution_ad = A*y_alpha_mantle_ad.y + D*y_beta_mantle_ad.y + E*y_gamma_mantle_ad.y

            #Combine solutions from all layers
            self.y_comb_solution = np.concatenate((y_inner_solution_ad, y_outer_solution_ad, y_mantle_solution_ad), axis=1)

            if verbose: print('Combined all solutions!')

            ###Love numbers calculation###
            self.h_n = self.y_comb_solution[0,-1]
            self.l_n = self.y_comb_solution[2,-1]
            self.k_n = self.y_comb_solution[4,-1] - 1
            self.delta_n = 1 + 2*self.h_n/self.n - (self.n+1)*self.k_n/self.n

            if verbose:
                print('Love numbers:')
                print(f'h_{n} = {self.h_n:.6f}')
                print(f'l_{n} = {self.l_n:.6f}')
                print(f'k_{n} = {self.k_n:.6f}')
                print(f'delta_{n} = {self.delta_n:.6f}')

        else:
            raise ValueError('We haven\'t implemented this setup yet.')


    #[PENDING] Put ouside integrate_internal_solutions_ad the setups
