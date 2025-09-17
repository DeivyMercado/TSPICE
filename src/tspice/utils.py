import numpy as np
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
