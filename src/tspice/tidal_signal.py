#Libraries
import numpy as np
import spiceypy as spy
import scipy.special as sps
import os
from tspice.utils import loc_func, convert_step_to_seconds
from tspice.kernels import download_kernels, write_meta_kernel

#Directories
file_dir = os.path.dirname(os.path.abspath(__file__)) #Current file location
data_dir = os.path.join(file_dir, 'data')	#Data files

#To track if kernels are loaded
kernels_loaded = False

def initialize(data_directory=None):

    """
    Download and load the necessary SPICE kernels.
    
    Input:
    - data_directory : [str, optional] Path to the directory where kernels should be stored and loaded from. If not provided, the default package data directory is used.
    """

    global kernels_loaded, data_dir
    
    #User can specify another data directory
    if data_directory is not None:
        data_dir = data_directory
    
    #Download kernels if not present
    download_kernels(data_dir)
    
    #Write the meta kernel and get its path
    meta_kernel_path = write_meta_kernel(data_dir)
    
    #Load the meta kernel
    spy.furnsh(meta_kernel_path)
    
    kernels_loaded = True
    print(f"TSPICE initialized successfully. Kernels loaded from: {meta_kernel_path}")

#This is the class that define the body
class Body():
	
	#Constructor of the class
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
			if not kernels_loaded:
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

	#Function to define array the times
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

	#Kn constants
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
	
	#Function to get the subpoint coordinates of an external body
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

	#Función del potencial de marea
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
	
	#Function to get the total tidal potential from multiple bodies
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