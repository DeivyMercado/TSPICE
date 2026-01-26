import numpy as np
import pandas as pd
import os
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from scipy.integrate import solve_ivp

#New density, mu and K
rho_new = 2178.9	#kg/m^3
mu_new = 266.6	#kbar
K_new = 520	#kbar
lam_new = K_new - (2/3)*mu_new	#kbar

#New velocities
Vp_new = 1e-3*np.sqrt((lam_new + 2*mu_new)*1e8/rho_new)	#m/s
Vs_new = 1e-3*np.sqrt(mu_new*1e8/rho_new)	#m/s
#print(Vp_new, Vs_new)

#Upload the data from the original paper by Dziewonski & Anderson (1981)
file_dir = os.path.dirname(os.path.abspath(__file__))
prem_paper = pd.read_excel(os.path.join(file_dir, 'PREM_amorin.xlsx'), sheet_name=0)
prem_paper.tail(10)

#Proposed modification: Top layer of the crust at 16 km depth
depth_mod = 16.0  #depth in km
r_top_crust = 6371.0 - depth_mod

#First split the PREM model into two parts: above and below the top of the crust
mask_bottom = prem_paper.loc[prem_paper['Radius']<=r_top_crust].index
prem_amorin_bottom = prem_paper.loc[mask_bottom]

#Change the units to SI
prem_amorin_bottom.Density = prem_amorin_bottom.Density*1e3		#Convert to kg/m^3
prem_amorin_bottom.K = prem_amorin_bottom.K*1e8		#Convert to Pa
prem_amorin_bottom.Mu = prem_amorin_bottom.Mu*1e8		#Convert to Pa
prem_amorin_bottom.Pressure = prem_amorin_bottom.Pressure*1e8	#Convert to Pa

#In both include the same radius to count for discontinuity
mask_top = prem_paper.loc[prem_paper['Radius']>=r_top_crust].index
prem_amorin_top = prem_paper.loc[mask_top]

#Transition layer
prem_amorin_top.loc[[85,86], 'Radius'] = r_top_crust

#Change the units to SI
prem_amorin_top.Density = prem_amorin_top.Density*1e3		#Convert to kg/m^3
prem_amorin_top.K = prem_amorin_top.K*1e8		#Convert to Pa
prem_amorin_top.Mu = prem_amorin_top.Mu*1e8		#Convert to Pa
prem_amorin_top.Pressure = prem_amorin_top.Pressure*1e8	#Convert to Pa

#Now change the values in the top part
prem_amorin_top.loc[86:, ['Density','Vp','Vs','K','Mu']] = [rho_new, Vp_new, Vs_new, K_new*1e8, mu_new*1e8]	
prem_amorin_top

#Erase previous discontinuities in the top layer
ind_discont_prev = prem_amorin_top[prem_amorin_top.duplicated(subset=['Radius'])].index.tolist()
prem_amorin_top.drop(ind_discont_prev[1:], inplace=True)

#Now concatenate both parts again
prem_amorin = pd.concat([prem_amorin_bottom, prem_amorin_top], ignore_index=True)

#Adding the Lamé parameter
prem_amorin['Lambda'] = prem_amorin['K'] - (2/3)*prem_amorin['Mu']	#Lambda in Pa

#Duplicated radius values (discontinuities)
ind_dupl_radius = prem_amorin[prem_amorin.duplicated(subset=['Radius'])].index.to_list()
n_dupl_radius = prem_amorin.duplicated(subset=['Radius']).sum()
r_dupl = prem_amorin.loc[ind_dupl_radius, 'Radius'].values
p_dupl = prem_amorin.loc[ind_dupl_radius, 'Pressure'].values

#Create interpolators for each layer
func = {}
for d in range(n_dupl_radius+1):
	if d == 0:
			ind_slice = slice(None, ind_dupl_radius[d])
	elif d == n_dupl_radius:
			ind_slice = slice(ind_dupl_radius[d-1], None)
	else:
			ind_slice = slice(ind_dupl_radius[d-1], ind_dupl_radius[d])
	
	#Extract data for the current layer and create interpolators
	r_d = prem_amorin['Radius'][ind_slice].values*1e3 #Radius in [m]
	lamb_d = prem_amorin['Lambda'][ind_slice].values
	mu_d = prem_amorin['Mu'][ind_slice].values
	p_d = prem_amorin['Pressure'][ind_slice].values
	rho_d = prem_amorin['Density'][ind_slice].values

	func[f'lamb_r_d{d}'] = CubicSpline(r_d, lamb_d, bc_type='natural')
	func[f'mu_r_d{d}'] = CubicSpline(r_d, mu_d, bc_type='natural')
	func[f'rho_r_d{d}'] = CubicSpline(r_d, rho_d, bc_type='natural')

	func[f'rho_P_d{d}'] = CubicSpline(p_d[::-1], rho_d[::-1], bc_type='natural')
	#func[f'rho_d{d}'] =interp1d(p_d, rho_d)

#Interpolator for density as a function of Pressure
def rho_P(P):
	idx = 7 - int(np.searchsorted(np.sort(p_dupl), P))
	return func[f'rho_P_d{idx}'](P)

#Interpolator for density as a function of radius
def rho_r(r):
	idx = int(np.searchsorted(r_dupl*1e3, r))
	return func[f'rho_r_d{idx}'](r)

#Interpolator for lambda as a function of radius
def lamb_r(r):
	idx = int(np.searchsorted(r_dupl*1e3, r))
	return func[f'lamb_r_d{idx}'](r)

#Interpolator for mu as a function of radius
def mu_r(r):
	idx = int(np.searchsorted(r_dupl*1e3, r))
	return func[f'mu_r_d{idx}'](r)

rho_P_interp = np.vectorize(rho_P)
rho_r_interp = np.vectorize(rho_r)
lamb_r_interp = np.vectorize(lamb_r)
mu_r_interp = np.vectorize(mu_r)

#Gravitational constant
G = 6.67430e-11	# m^3 kg^-1 s^-2

#Mass enclosed function
def mr(r, rho_interp):
  m = 4*np.pi*r**2*rho_interp(r)
  return m

#Equation of Hydrostatic Equilibrium and Mass Conservation
def edo_hidrostatico(r, Y, rho_interp):
  p, Mr = Y
  dpdr = -G*Mr*rho_interp(p)/r**2
  dMrdr = 4*np.pi*r**2*rho_interp(p)
  return [dpdr, dMrdr]

#Integration to get the total mass of the Earth in this model
M, dM = quad(mr, 0, 6371e3, args=(rho_r_interp,), epsrel=1e-3)
Mearth = M

#Radius array for the integration results
rs_integration = np.linspace(6e3, 6371e3, 100000)
rs_integration_mass = rs_integration[::-1]

#Integration of the EDOs
Y0 = [0, Mearth]
solution = solve_ivp(edo_hidrostatico, (rs_integration_mass[0], rs_integration_mass[-1]), Y0, t_eval=rs_integration_mass, args=(rho_P_interp,), rtol=1e-6, atol=1e-9)
ps = solution.y[0]
Mrs = solution.y[1]

#Interpolation of gravity excluding the innermost part and forcing to zero in the center
Mrs_c = Mrs[::-1] #from center to surface
Mrs_redu = Mrs_c[rs_integration>1000e3]
rs_redu = rs_integration[rs_integration>1000e3]
gs_redu = G*Mrs_redu/(rs_redu**2)
x = np.concatenate([[0.0], rs_redu.flatten()])
y = np.concatenate([[0.0], gs_redu.flatten()])
g_r_interp = CubicSpline(x, y)

#Mass interpolator
M_r_interp = CubicSpline(rs_integration, Mrs_c)
