import numpy as np

###Adimensional differential equations

#TAKE CARE: In the definition of the differential equations, uses 'y' first and 'r' second in odeint, and 'r' first and 'y' second in solve_ivp

#Solid layers:

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

#Fluid layers:

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
