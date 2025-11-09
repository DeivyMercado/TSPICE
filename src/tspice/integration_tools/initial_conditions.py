import numpy as np

### Initial conditions

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

    #Auxiliar variable
    J = mu**(-1)
    gamma = (4/3)*np.pi*rho
    
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