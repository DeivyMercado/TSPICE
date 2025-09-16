import numpy as np
import re

#To get coordinates for the Vtid/g expression
def loc_func(loc, a_main):

    '''This function converts the geographic coordinates of a station (lon, lat, depth) to spherical coordinates (phi, theta, a)
    '''

    #Coordinates of the station
    lon_sta, lat_sta, depth_sta = loc['lon'], loc['lat'], loc['depth']
    colat_sta = 90 - lat_sta    #Colatitud
    theta_sta, phi_sta = np.deg2rad(colat_sta), np.deg2rad(lon_sta)

    #[PENDIENT] Revisar esto de la distancia entre el ellipsoide y la estación
    #Distance from the COM to the station
    a_sta = a_main - depth_sta

    return phi_sta, theta_sta, a_sta

#Matching regular expressions in 'step' to get the step in seconds
def convert_step_to_seconds(step):
    '''This function converts the time step from regular expressions like "15s", "30m", "1h", "2d", etc. to seconds.

    Input:
    - step: Time step as #s, #m, #h, #d, etc.

    Output:
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
