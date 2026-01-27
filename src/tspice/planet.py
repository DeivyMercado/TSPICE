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
Planetary Models Module.

This module provides planetary interior models for use in tidal response calculations.
Currently implements the PREM (Preliminary Reference Earth Model) with modifications
from Amorin & Gudkova (2025).

Classes
-------
Earth
    Earth interior model based on PREM with modified crust properties.

Examples
--------
Create an Earth model and access its properties:

>>> from tspice import Earth
>>> earth = Earth()
>>> profile = earth.planet_profile
>>> print(f"Density at center: {profile['rho'](0):.2f} kg/m^3")

Use Earth model for Love number calculations (integration):

>>> import tspice as tsp
>>> from tspice.planet import Earth
>>> 
>>> # 1. Define responding body and characteristic scales
>>> earth_interior = tsp.BodyResponse('Earth')
>>> earth_interior.scale_constants()
>>> 
>>> # 2. Load planet profile (Modified PREM)
>>> earth_model = Earth()
>>> planet_profile = earth_model.planet_profile
>>> 
>>> # 3. Define internal layers (Fluid Outer Core sandwiched between Solid Inner Core and Mantle)
>>> layers_list = [
...     dict(name='Outer Core', type='fluid', r0=1221500.0, rf=3480000.0),
...     dict(name='Inner Core', type='solid', r0=0, rf=1221500.0),
...     dict(name='Mantle + crust', type='solid', r0=3480000.0, rf=earth_interior.L)
... ]
>>> 
>>> # 4. Set integration parameters (e.g., M2 tide, n=2)
>>> f_day = 1.93502
>>> earth_interior.set_integration_parameters_ad(
...     n=2, f_days=f_day,
...     layers_list=layers_list,
...     planet_profile=planet_profile,
...     r0_ini=6e3
... )
>>> 
>>> # 5. Integrate coupled ODEs
>>> earth_interior.integrate_internal_solutions_ad(verbose=True)
>>> print(f"k_2 = {earth_interior.k_n}, h_2 = {earth_interior.h_n}, l_2 = {earth_interior.l_n}")
"""

import numpy as np
import pandas as pd
import os
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from scipy.integrate import solve_ivp
from tspice import get_data

class Earth:
    def __init__(self):
        """
        Initialize the Earth interior model.

        This method loads the PREM model data, applies the crustal modifications
        suggested by Amorin & Gudkova (2025), and calculates the gravity profile
        and mass distribution required for tidal response integration.
        """
        # New density, mu and K parameters
        self.rho_new = 2178.9  # kg/m^3
        self.mu_new = 266.6    # kbar
        self.K_new = 520       # kbar
        self.lam_new = self.K_new - (2/3) * self.mu_new  # kbar

        # New velocities
        self.Vp_new = 1e-3 * np.sqrt((self.lam_new + 2 * self.mu_new) * 1e8 / self.rho_new)  # m/s
        self.Vs_new = 1e-3 * np.sqrt(self.mu_new * 1e8 / self.rho_new)  # m/s

        # Constants
        self.G = 6.67430e-11  # m^3 kg^-1 s^-2
        self.Re = 6371e3      # Radius of Earth in m

        # Initialize data processing and interpolators
        self._initialize_model()

    def _initialize_model(self):
        # Load PREM data
        prem_file = get_data("PREM_amorin.xlsx")
        prem_paper = pd.read_excel(prem_file, sheet_name=0)

        # Proposed modification: Top layer of the crust at 16 km depth
        depth_mod = 16.0  # depth in km
        r_top_crust = 6371.0 - depth_mod

        # Split PREM model
        mask_bottom = prem_paper.loc[prem_paper['Radius'] <= r_top_crust].index
        prem_amorin_bottom = prem_paper.loc[mask_bottom].copy()

        # Change units to SI for bottom part
        prem_amorin_bottom.loc[:, 'Density'] *= 1e3
        prem_amorin_bottom.loc[:, 'K'] *= 1e8
        prem_amorin_bottom.loc[:, 'Mu'] *= 1e8
        prem_amorin_bottom.loc[:, 'Pressure'] *= 1e8

        mask_top = prem_paper.loc[prem_paper['Radius'] >= r_top_crust].index
        prem_amorin_top = prem_paper.loc[mask_top].copy()

        # Transition layer
        prem_amorin_top.loc[[85, 86], 'Radius'] = r_top_crust

        # Change units to SI for top part (before replacing modification values)
        prem_amorin_top.loc[:, 'Density'] *= 1e3
        prem_amorin_top.loc[:, 'K'] *= 1e8
        prem_amorin_top.loc[:, 'Mu'] *= 1e8
        prem_amorin_top.loc[:, 'Pressure'] *= 1e8

        # Cast columns to float before assignment to avoid TypeError
        cols_to_update = ['Density', 'Vp', 'Vs', 'K', 'Mu']
        prem_amorin_top[cols_to_update] = prem_amorin_top[cols_to_update].astype(float)

        # Apply modification values
        prem_amorin_top.loc[86:, cols_to_update] = [
            self.rho_new, self.Vp_new, self.Vs_new, self.K_new * 1e8, self.mu_new * 1e8
        ]

        # Erase previous discontinuities in top layer
        ind_discont_prev = prem_amorin_top[prem_amorin_top.duplicated(subset=['Radius'])].index.tolist()
        if len(ind_discont_prev) > 1:
             prem_amorin_top.drop(ind_discont_prev[1:], inplace=True)

        # Concatenate
        self.prem_amorin = pd.concat([prem_amorin_bottom, prem_amorin_top], ignore_index=True)

        # Add Lambda parameter
        self.prem_amorin['Lambda'] = self.prem_amorin['K'] - (2/3) * self.prem_amorin['Mu']

        # Determine discontinuities
        ind_dupl_radius = self.prem_amorin[self.prem_amorin.duplicated(subset=['Radius'])].index.to_list()
        n_dupl_radius = len(ind_dupl_radius)
        
        self.r_dupl = self.prem_amorin.loc[ind_dupl_radius, 'Radius'].values
        self.p_dupl = self.prem_amorin.loc[ind_dupl_radius, 'Pressure'].values

        # Creating interpolators
        self._func = {}
        for d in range(n_dupl_radius + 1):
            if d == 0:
                ind_slice = slice(None, ind_dupl_radius[d])
            elif d == n_dupl_radius:
                ind_slice = slice(ind_dupl_radius[d-1], None)
            else:
                ind_slice = slice(ind_dupl_radius[d-1], ind_dupl_radius[d])

            r_d = self.prem_amorin['Radius'][ind_slice].values * 1e3
            lamb_d = self.prem_amorin['Lambda'][ind_slice].values
            mu_d = self.prem_amorin['Mu'][ind_slice].values
            p_d = self.prem_amorin['Pressure'][ind_slice].values
            rho_d = self.prem_amorin['Density'][ind_slice].values

            self._func[f'lamb_r_d{d}'] = CubicSpline(r_d, lamb_d, bc_type='natural')
            self._func[f'mu_r_d{d}'] = CubicSpline(r_d, mu_d, bc_type='natural')
            self._func[f'rho_r_d{d}'] = CubicSpline(r_d, rho_d, bc_type='natural')
            self._func[f'rho_P_d{d}'] = CubicSpline(p_d[::-1], rho_d[::-1], bc_type='natural')

        # Create vectorized interpolators
        self.rho_P_interp = np.vectorize(self._rho_P)
        self.rho_r_interp = np.vectorize(self._rho_r)
        self.lamb_r_interp = np.vectorize(self._lamb_r)
        self.mu_r_interp = np.vectorize(self._mu_r)

        # Calculate Gravity Profile
        self._calculate_gravity()

    # Interpolator Methods
    def _rho_P(self, P):
        idx = 7 - int(np.searchsorted(np.sort(self.p_dupl), P))
        return self._func[f'rho_P_d{idx}'](P)

    def _rho_r(self, r):
        idx = int(np.searchsorted(self.r_dupl * 1e3, r))
        return self._func[f'rho_r_d{idx}'](r)

    def _lamb_r(self, r):
        idx = int(np.searchsorted(self.r_dupl * 1e3, r))
        return self._func[f'lamb_r_d{idx}'](r)

    def _mu_r(self, r):
        idx = int(np.searchsorted(self.r_dupl * 1e3, r))
        return self._func[f'mu_r_d{idx}'](r)

    def _calculate_gravity(self):
        # Mass enclosed function
        def mr(r, rho_interp):
            return 4 * np.pi * r**2 * rho_interp(r)

        # Integration for total mass
        M, dM = quad(mr, 0, self.Re, args=(self.rho_r_interp,), epsrel=1e-3)
        self.Mearth = M

        # Equation of Hydrostatic Equilibrium
        def edo_hidrostatico(r, Y, rho_interp):
            p, Mr = Y
            dpdr = -self.G * Mr * rho_interp(p) / r**2
            dMrdr = 4 * np.pi * r**2 * rho_interp(p)
            return [dpdr, dMrdr]

        rs_integration = np.linspace(6e3, self.Re, 100000)
        rs_integration_mass = rs_integration[::-1]

        Y0 = [0, self.Mearth]
        solution = solve_ivp(edo_hidrostatico, (rs_integration_mass[0], rs_integration_mass[-1]), Y0, 
                             t_eval=rs_integration_mass, args=(self.rho_P_interp,), rtol=1e-6, atol=1e-9)
        
        Mrs = solution.y[1]
        Mrs_c = Mrs[::-1] # from center to surface

        # Interpolation of gravity
        mask = rs_integration > 1000e3
        Mrs_redu = Mrs_c[mask]
        rs_redu = rs_integration[mask]
        gs_redu = self.G * Mrs_redu / (rs_redu**2)

        x = np.concatenate([[0.0], rs_redu.flatten()])
        y = np.concatenate([[0.0], gs_redu.flatten()])
        
        self.g_r_interp = CubicSpline(x, y)
        self.M_r_interp = CubicSpline(rs_integration, Mrs_c)

    @property
    def planet_profile(self):
        """
        Returns the planet profile dictionary compatible with BodyResponse.

        The profile contains vectorized interpolators for physical properties
        as a function of radius `r` (in meters).

        Returns
        -------
        dict
            Dictionary containing the following keys (values are callable interpolators):
            - 'rho' : Density [kg/m^3]
            - 'lamb': First Lamé parameter [Pa]
            - 'mu'  : Shear modulus (Second Lamé parameter) [Pa]
            - 'g'   : Gravitational acceleration [m/s^2]
            - 'dimensionless': Boolean flag (False for this dimensional model)
        
        Examples
        --------
        >>> earth = Earth()
        >>> profile = earth.planet_profile
        >>> rho_surface = profile['rho'](6371000)
        >>> g_surface = profile['g'](6371000)
        """
        return {
            'rho': self.rho_r_interp,
            'lamb': self.lamb_r_interp,
            'mu': self.mu_r_interp,
            'g': self.g_r_interp,
            'dimensionless': False
        }
