##################################################################
#                                                                #
# tSPICE: Tidal Signal with Python and SPICE                     #
#                                                                #
##################################################################
# License: GNU Affero General Public License v3 (AGPL-3.0)       #
# Authors: Jorge I. Zuluaga, Juanita A. Agudelo                  #
# Contact: jorge.zuluaga@udea.edu.co                             #
##################################################################

import numpy as np
import spiceypy as spy
from datetime import datetime
import re

#To get the triaxial ellipsoid radius
def ellipsoid_radius(phi, theta, radii):
    
    """
    Calculate the radial distance to the surface of the triaxial ellipsoid.
    
    Inputs:
    - phi: [float or array-like] Longitude in radians (0 at prime meridian)
	- theta: [float or array-like] Colatitude in radians (0 at North pole, π at South pole)
    - a, b, c: [float] Semi-axes of the triaxial ellipsoid (typically a ≥ b ≥ c)

    Outputs:
    - r : [float or array-like] Radial distance from center to surface at (phi, theta)
    """

    a, b, c = radii
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    
    #Formula for triaxial ellipsoid radius
    r = 1/np.sqrt((sin_theta * cos_phi / a)**2
                          + (sin_theta * sin_phi / b)**2
                          + (cos_theta / c)**2)
    
    return r

#To get coordinates for the Vtid/g expression
def loc_func(loc, radii):

    '''
    This function converts the geographic coordinates of a station (lon, lat, depth) to spherical coordinates (phi, theta, a).

    Inputs:
    - loc: [deg, deg, km] Dictionary with the geographic coordinates of the station. For instance: dict(lon = -70.0, lat = 40.0, depth = 10.0)
    - radii: [km] Array with the semi-axes of the triaxial ellipsoid (a, b, c).

    Outputs:
    - phi_sta: [rad] Longitude in radians.
    - theta_sta: [rad] Colatitude in radians.
    - a_sta: [km] Ellipsoidal radius at the station coordinates.
    '''

    #Coordinates of the station
    lon_sta, lat_sta, depth_sta = loc['lon'], loc['lat'], loc['depth']
    colat_sta = 90 - lat_sta    #Colatitude
    phi_sta, theta_sta = np.deg2rad(lon_sta), np.deg2rad(colat_sta)

    #Ellipsoidal radius at the station coordinates
    a_ellip = ellipsoid_radius(phi_sta, theta_sta, radii)
    
    #Distance from the COM to the station
    a_sta = a_ellip - depth_sta     #Assuming the depth is respect to the ellipsoidal of reference

    return phi_sta, theta_sta, a_sta

#Matching regular expressions in 'step' to get the step in seconds
def convert_step_to_seconds(step):

    '''
    This function converts the time step from regular expressions like "15s", "30m", "1h", "2d", etc. to seconds.

    Inputs:
    - step: Time step as #s, #m, #h, #d, etc.

    Outputs:
    - step_s: Step in seconds.
    '''

    #Match number and unit
    matched = re.match(r"(\d+)([smhd])", step)
    if not matched:
        raise ValueError("Invalid step format. Use formats like '1h', '30s', etc.")
    value, unit = int(matched.group(1)), matched.group(2)

    #Conversion dictionary
    unit_conversion = {'s': 1,    	#Segundos
                    'm': 60,   	#Minutos
                    'h': 3600, 	#Horas
                    'd': 86400}	#Días

    #Conversion to seconds
    step_s = value*unit_conversion[unit]

    return step_s

def et_to_utc_string(et, prec=2, format='ISOC', datetime_format=True):
    
    '''
    This function converts Ephemeris Time (ET) to a UTC string format "YYYY-MM-DD HR:MN:SC".

    Inputs:
    - et: [float] Ephemeris Time.

    Outputs:
    - utc_string: [str] UTC time in the format "YYYY-MM-DD HR:MN:SC".
    '''

    #Convert ET to UTC with spiceypy
    utc_times = spy.et2utc(et, format_str=format, prec=prec)	#UTC strings

    #Convert to datetime objects if needed (useful for plotting)
    if datetime_format:
        utc_times = [datetime.fromisoformat(utc) for utc in utc_times] #Datetime objects

    return utc_times

import numpy as np



def Y0_TakeuchiSaito1972_ad(r0, params_ad):
    '''
    Set of three independent solutions for y_i near the center of the Earth in adimensional form. From the power series expansion in Crossley (1975), using the definitions in Takeuchi & Saito (1972).
        
    Input:
    - r0 : Small radius near the center of the Earth.
    - params : Dictionary containing: n, lam, mu, rho, g.

    where:
    - n: degree of the tidal potential
    - lam, mu: Lamé parameters
    - rho: density at equilibrium
    - g: gravitational acceleration at equilibrium
    
    Output:
    - Y0 : array-like Vector of dependent variables [y1, y2, y3, y4, y5, y6] at r=r0.
    '''
    
    #Degree of the tidal potential
    n = params_ad["n"]
    
    #Planetary profiles at r=r0
    lam = params_ad["lam"](0)
    mu = params_ad["mu"](0)
    rho = params_ad["rho"](0)

	#First independent solution
    y11 = r0**(n-1)
    y21 = (lam*(n+1) + 2*mu*(n-1))*r0**(n-2)
    y31 = 0
    y41 = mu*r0**(n-2)
    y51 = 0
    y61 = -4*np.pi*rho*r0**(n-1)
    
	#Second independent solution
    y12 = 0
    y22 = -lam*n*(n+1)*r0**(n-2)
    y32 = r0**(n-1)
    y42 = mu*(n-2)*r0**(n-2)
    y52 = 0
    y62 = 0

	#Third independent solution
    y13 = 0
    y23 = 0
    y33 = 0
    y43 = 0
    y53 = r0**n
    y63 = (2*n+1)*r0**(n-1)

    return np.array([[y11, y21, y31, y41, y51, y61],
                    [y12, y22, y32, y42, y52, y62],
					[y13, y23, y33, y43, y53, y63]]).transpose()

def Y0_XuSun2003_ad(r0, params_ad):
    '''
    Set of three independent solutions for y_i near the center of the Earth in adimensional form. From the power series expansion in Crossley (1975) and Xu & Sun (2003).
        
    Input:
    - r0 : Small radius near the center of the Earth.
    - params : Dictionary containing: n, lam, mu, rho, g.

    where:
    - n: degree of the tidal potential
    - lam, mu: Lamé parameters
    - rho: density at equilibrium
    - g: gravitational acceleration at equilibrium
    
    Output:
    - Y0 : array-like Vector of dependent variables [y1, y2, y3, y4, y5, y6] at r=r0.
    '''
    
    #Degree of the tidal potential
    n = params_ad["n"]
    
    #Planetary profiles at r=r0
    lam = params_ad["lam"](0)
    mu = params_ad["mu"](0)
    rho = params_ad["rho"](0)
    
	#First independent solution
    y11 = r0**(n-1)
    y21 = (lam*(n+1) + 2*mu*(n-1))*r0**(n-2)
    y31 = 0
    y41 = mu*r0**(n-2)
    y51 = 0
    y61 = 4*np.pi*rho*r0**(n-1)
    
	#Second independent solution
    y12 = 0
    y22 = -lam*n*(n+1)*r0**(n-2)
    y32 = r0**(n-1)
    y42 = mu*(n-2)*r0**(n-2)	#This was an error in Xu & Sun (2003)
    y52 = 0
    y62 = 0

	#Third independent solution
    y13 = 0
    y23 = 0
    y33 = 0
    y43 = 0
    y53 = r0**n
    y63 = n*r0**(n-1)

    return np.array([[y11, y21, y31, y41, y51, y61],
                    [y12, y22, y32, y42, y52, y62],
					[y13, y23, y33, y43, y53, y63]]).transpose()

def Y0_AmorinGudkova2024_ad(r0, params_ad):
    '''
    Set of three independent solutions for y_i near the center of the Earth in adimensional form. From the power series expansion in Amorin & Gudkova (2024). Compared to Xu & Sun (2003), includes higher order terms in y_1 and y_3.
        
    Input:
    - r0 : Small radius near the center of the Earth.
    - params : Dictionary containing: n, lam, mu, rho, g.

    where:
    - n: degree of the tidal potential
    - omega: frequency parameter
    - lam, mu: Lamé parameters
    - rho: density at equilibrium
    - g: gravitational acceleration at equilibrium
    
    Output:
    - Y0 : array-like Vector of dependent variables [y1, y2, y3, y4, y5, y6] at r=r0.
    '''
    
    #Degree of the tidal potential
    n = params_ad["n"]
    omega = params_ad["omega"]
    
    #Planetary profiles at r=r0
    lam = params_ad["lam"](0)
    mu = params_ad["mu"](0)
    rho = params_ad["rho"](0)

    #Auxiliar variable
    gamma = (4/3)*np.pi*rho
    J = gamma/(n*gamma - omega**2)
    
	#First independent solution
    y11 = n*J*r0**(n-1)
    y21 = 2*n*(n-1)*r0**(n-2)
    y31 = J*r0**(n-1)
    y41 = 2*(n-1)*r0**(n-2)
    y51 = gamma*r0**n
    y61 = gamma*n*(1-3*J)*r0**(n-1)
    
	#Second independent solution
    y12 = 0
    y22 = -n*(n+1)*lam*r0**n
    y32 = r0**(n+1)
    y42 = n*mu*r0**n
    y52 = 0
    y62 = 0

	#Third independent solution
    y13 = r0**(n+1)
    y23 = (n*(lam + 2*mu) + 3*lam + 2*mu)*r0**n
    y33 = 0
    y43 = mu*r0**n
    y53 = 0
    y63 = -3*gamma*r0**(n+1)

    return np.array([[y11, y21, y31, y41, y51, y61],
                    [y12, y22, y32, y42, y52, y62],
					[y13, y23, y33, y43, y53, y63]]).transpose()



#TAKE CARE: In the definition of the differential equations, uses 'y' first and 'r' second in odeint, and 'r' first and 'y' second in solve_ivp



def dydr_solid_TakeuchiSaito1972_ad(r, y, params_adim):
		
		"""
		Differential equations for the adimensional internal solution of an Earth-like planet interior in the solid layers, using the convention of Takeuchi & Saito (1972).

		Input:
		- r: radius.
		- y: Vector of dependent variables [y1, y2, y3, y4, y5, y6].
		- params_adim: Dictionary containing the adiminensional parameters:

		where:
		- n: degree of the tidal potential
		- omega: angular frequency of the tidal potential
		- lam, mu: Lamé parameters
		- rho: density at equilibrium
		- g: gravitational acceleration at equilibrium
		
		Output:
        - dy/dr: Vector of derivatives of the dependent variables [dy1/dr, dy2/dr, dy3/dr, dy4/dr, dy5/dr, dy6/dr].
		"""

		#Unpack variables
		y1, y2, y3, y4, y5, y6 = y

		#Degree of the tidal potential
		n = params_adim["n"]
		n1 = n*(n + 1)

		#Angular frequency
		omega = params_adim["omega"]

		#Planetary profiles
		lam = params_adim["lam"](r)
		mu = params_adim["mu"](r)
		rho = params_adim["rho"](r)
		g = params_adim["g"](r)
		beta = lam + 2*mu
		gamma = mu*(3*lam + 2*mu)/beta
		
        #Differential equations
		dy1 = -2*lam*y1/(beta*r) + y2/beta + n1*lam*y3/(beta*r)
		dy2 = -(rho*omega**2*r**2 + 4*g*rho*r - 4*gamma)*y1/r**2 - 4*mu*y2/(beta*r) + (n1*g*rho*r - 2*n1*gamma)*y3/r**2 + n1*y4/r + rho*y6
		dy3 = -y1/r + y3/r + y4/mu
		dy4 = (g*rho*r - 2*gamma)*y1/r**2 - lam*y2/(beta*r) + (-rho*omega**2*r**2 + (2*mu/beta)*((2*n1 - 1)*lam + 2*(n1 - 1)*mu))*y3/r**2 - 3*y4/r - rho*y5/r
		dy5 = 4*np.pi*rho*y1 - (n+1)*y5/r + y6
		dy6 = 4*np.pi*rho*(n+1)*y1/r - 4*np.pi*rho*n1*y3/r + (n-1)*y6/r

		return np.array([dy1, dy2, dy3, dy4, dy5, dy6])

def dydr_solid_XuSun2003_ad(r, y, params_adim):
		
		"""
		Differential equations for the adimensional internal solution of an Earth-like planet interior in the solid layers, using the convention of Xu & Sun (2003).

		Input:
		- r: radius.
		- y: Vector of dependent variables [y1, y2, y3, y4, y5, y6].
		- params_adim: Dictionary containing the adiminensional parameters:

		where:
		- n: degree of the tidal potential
		- omega: angular frequency of the tidal potential
		- lam, mu: Lamé parameters
		- rho: density at equilibrium
		- g: gravitational acceleration at equilibrium
		
		Output:
        - dy/dr: Vector of derivatives of the dependent variables [dy1/dr, dy2/dr, dy3/dr, dy4/dr, dy5/dr, dy6/dr].
		"""

		#Unpack variables
		y1, y2, y3, y4, y5, y6 = y

		#Degree of the tidal potential
		n = params_adim["n"]
		n1 = n*(n + 1)

		#Angular frequency
		omega = params_adim["omega"]

		#Planetary profiles
		lam = params_adim["lam"](r)
		mu = params_adim["mu"](r)
		rho = params_adim["rho"](r)
		g = params_adim["g"](r)
		beta = lam + 2*mu
		gamma = mu*(3*lam + 2*mu)/beta

        #Differential equations
		dy1 = -2*lam*y1/(beta*r) + y2/beta + n1*lam*y3/(beta*r)
		dy2 = -(rho*omega**2*r**2 + 4*g*rho*r - 4*gamma)*y1/r**2 - 4*mu*y2/(beta*r) + (n1*g*rho*r - 2*n1*gamma)*y3/r**2 + n1*y4/r + rho*y6
		dy3 = -y1/r + y3/r + y4/mu
		dy4 = (g*rho*r - 2*gamma)*y1/r**2 - lam*y2/(beta*r) + (-rho*omega**2*r**2 + (2*mu/beta)*((2*n1 - 1)*lam + 2*(n1 - 1)*mu))*y3/r**2 - 3*y4/r - rho*y5/r
		dy5 = -4*np.pi*rho*y1 + y6
		dy6 = 4*np.pi*rho*n1*y3/r + n1*y5/r**2 - 2*y6/r

		return np.array([dy1, dy2, dy3, dy4, dy5, dy6])

def dydr_solid_AmorinGudkova2024_ad(r, y, params_adim):
		
		"""
		Differential equations for the adimensional internal solution of an Earth-like planet interior in the solid layers, using the convention of Amorin & Gudkova (2024).

		Input:
		- r: radius.
		- y: Vector of dependent variables [y1, y2, y3, y4, y5, y6].
		- params_adim: Dictionary containing the adiminensional parameters:

		where:
		- n: degree of the tidal potential
		- omega: angular frequency of the tidal potential
		- lam, mu: Lamé parameters
		- rho: density at equilibrium
		- g: gravitational acceleration at equilibrium
		
		Output:
        - dy/dr: Vector of derivatives of the dependent variables [dy1/dr, dy2/dr, dy3/dr, dy4/dr, dy5/dr, dy6/dr].
		"""

		#Unpack variables
		y1, y2, y3, y4, y5, y6 = y

		#Degree of the tidal potential
		n = params_adim["n"]
		n1 = n*(n + 1)

		#Angular frequency
		omega = params_adim["omega"]

		#Planetary profiles
		lam = params_adim["lam"](r)
		mu = params_adim["mu"](r)
		rho = params_adim["rho"](r)
		g = params_adim["g"](r)
		beta = lam + 2*mu
		gamma = mu*(3*lam + 2*mu)/beta

        #Differential equations
		dy1 = -2*lam*y1/(beta*r) + y2/beta + n1*lam*y3/(beta*r)
		dy2 = -(rho*(omega**2)*(r**2) + 4*g*rho*r - 4*gamma)*y1/(r**2) - 4*mu*y2/(beta*r) + (n1*g*rho*r - 2*n1*gamma)*y3/(r**2) + n1*y4/r - rho*y6
		dy3 = -y1/r + y3/r + y4/mu
		dy4 = (g*rho*r - 2*gamma)*y1/(r**2) - lam*y2/(beta*r) + (-rho*omega**2*r**2 + (2*mu/beta)*((2*n1 - 1)*lam + 2*(n1 - 1)*mu))*y3/(r**2) - 3*y4/r - rho*y5/r	#The original paper had a typo here with a plus sign before rho*omega**2*r**2
		dy5 = 4*np.pi*rho*y1 + y6
		dy6 = -4*np.pi*rho*n1*y3/r + n1*y5/r**2 - 2*y6/r

		return np.array([dy1, dy2, dy3, dy4, dy5, dy6])



def dydr_fluid_TakeuchiSaito1972_ad(r, y, params_ad):
		
		"""
		Differential equations for the internal solution of an Earth-like planet interior in the fluid layers.

		Input:
		- r: Radius.
		- y: Vector of dependent variables [y1, y2, y3, y4, y5, y6].
		- params_ad: Dictionary containing:

		where:
		- n: degree of the tidal potential
		- omega: angular frequency of the tidal potential
		- lam: Lamé parameter
		- rho: density at equilibrium
		- g: gravitational acceleration at equilibrium
		"""

		#Unpack variables
		y1, y2, y5, y6 = y

		#Degree of the tidal potential
		n = params_ad["n"]
		n1 = n*(n + 1)

		#Angular frequency
		omega = params_ad["omega"]

		#Planetary profiles
		lam = params_ad["lam"](r)
		rho = params_ad["rho"](r)
		g = params_ad["g"](r)

		#Fixed values
		y3 = (g*y1 - y2/rho + y5)/(omega**2*r)
		y4 = 0

		#Liquid region (μ = 0)
		dy1 = -2*y1/r + y2/lam + n1*y3/r
		dy2 = -(rho*omega**2 + 4*rho*g/r)*y1 + n1*rho*g*y3/r + rho*y6
		dy5 = 4*np.pi*rho*y1 - (n+1)*y5/r + y6
		dy6 = 4*np.pi*rho*(n+1)*y1/r - 4*np.pi*rho*n1*y3/r + (n-1)*y6/r 

		return np.array([dy1, dy2, dy5, dy6])


def dydr_fluid_XuSun2003_ad(r, y, params_ad):
		
		"""
		Differential equations for the internal solution of an Earth-like planet interior in the fluid layers.

		Input:
		- r: Radius.
		- y: Vector of dependent variables [y1, y2, y3, y4, y5, y6].
		- params_ad: Dictionary containing:

		where:
		- n: degree of the tidal potential
		- omega: angular frequency of the tidal potential
		- lam: Lamé parameter
		- rho: density at equilibrium
		- g: gravitational acceleration at equilibrium
		"""

		#Unpack variables
		y1, y2, y5, y6 = y

		#Degree of the tidal potential
		n = params_ad["n"]
		n1 = n*(n + 1)

		#Angular frequency
		omega = params_ad["omega"]

		#Planetary profiles
		lam = params_ad["lam"](r)
		rho = params_ad["rho"](r)
		g = params_ad["g"](r)

		#Fixed values
		y3 = (g*y1 - y2/rho + y5)/(omega**2*r)
		y4 = 0

		#Liquid region (μ = 0)
		dy1 = -2*y1/r + y2/lam + n1*y3/r
		dy2 = -(rho*omega**2 + 4*rho*g/r)*y1 + n1*rho*g*y3/r + rho*y6
		dy5 = -4*np.pi*rho*y1 + y6
		dy6 = 4*np.pi*rho*n1*y3/r + n1*y5/r**2 - 2*y6/r

		return np.array([dy1, dy2, dy5, dy6])

def dydr_fluid_AmorinGudkova2024_ad(r, y, params_ad):
		
		"""
		Differential equations for the internal solution of an Earth-like planet interior in the fluid layers.
		Key difference with the Xu & Sun (2003) approach: sign in the equation for dy5/dr.

		Input:
		- r: Radius.
		- y: Vector of dependent variables [y1, y2, y3, y4, y5, y6].
		- params_ad: Dictionary containing:

		where:
		- n: degree of the tidal potential
		- omega: angular frequency of the tidal potential
		- lam: Lamé parameter
		- rho: density at equilibrium
		- g: gravitational acceleration at equilibrium
		"""

		#Unpack variables
		y1, y2, y5, y6 = y

		#Degree of the tidal potential
		n = params_ad["n"]
		n1 = n*(n + 1)

		#Angular frequency
		omega = params_ad["omega"]

		#Planetary profiles
		lam = params_ad["lam"](r)
		rho = params_ad["rho"](r)
		g = params_ad["g"](r)

		#Fixed values
		y3 = (g*y1 - y2/rho + y5)/(omega**2*r)
		y4 = 0

		#Liquid region (μ = 0)
		dy1 = -2*y1/r + y2/lam + n1*y3/r
		dy2 = -(rho*omega**2 + 4*rho*g/r)*y1 + n1*rho*g*y3/r + rho*y6
		dy5 = 4*np.pi*rho*y1 + y6
		dy6 = -4*np.pi*rho*n1*y3/r + n1*y5/r**2 - 2*y6/r

		return np.array([dy1, dy2, dy5, dy6])

def dzdr_fluid_AmorinGudkova2024_ad(r, z, params_ad):
		
		"""
		Differential equations for the internal solution of an Earth-like planet interior in the fluid layers.
		Key difference with the Xu & Sun (2003) approach: sign in the equation for dy5/dr.

		Input:
		- r: Radius.
		- z: Vector of dependent variables [z5, z7].
		- params_ad: Dictionary containing:

		where:
		- n: degree of the tidal potential
		- rho: density at equilibrium
		- g: gravitational acceleration at equilibrium
		"""

		#Unpack variables
		z5, z7 = z

		#Degree of the tidal potential
		n = params_ad["n"]

		#Planetary profiles
		rho = params_ad["rho"](r)
		g = params_ad["g"](r)

		#Simplified equations in the Liquid region (μ = 0)
		dz5 = (4*np.pi*rho/g - (n+1)/r)*z5 + z7
		dz7 = 2*(n-1)*4*np.pi*rho*z5/(r*g) + ((n-1)/r - 4*np.pi*rho/g)*z7

		return np.array([dz5, dz7])
