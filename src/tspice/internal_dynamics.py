#Libraries
import numpy as np
from scipy.integrate import solve_ivp

#Import TSPICE
from tspice import Body
from tspice.integration_tools.initial_conditions import *
from tspice.integration_tools.differential_equations import *

#This class will contain the body information to do the integration of equations of motion
class InteriorModel(Body):

    #Constructor of the class
    def __init__(self, name):

        #The class inherit the attributes from Body
        super().__init__(name)
        #GM_main, a_ellips, a_main, f_main, g_ref

    #To calculate the characteristic scales for adimensionalization
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
        self.L = int(self.a_main)*1e3	#Length in [m] 

        self.M = self.GM_main*1e9/self.G	#Mass in [kg]
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


    #To read the layers and define the steps for each layer
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
    
    #Integration parameters
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

    #Initial conditions for the integration
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

    #Integrate the yis functions
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
        