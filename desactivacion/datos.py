import numpy as np

# diseño del intercambiador
Dint = 5*2.54/100 # metros, diametro interno de los tubos
Ntub = 2400 # nuemro de tubos
L = 2.5 # metros, longitud de los tubos



# catalizador
ro_l = 1700. #densidad del lecho kg/m3_lecho
ro_p = 2000. #densidad del catalizador kg/m3_cat
dp = 3.e-3 #diametro medio del Catalizdor (m)
e = 1-ro_l/ro_p #porosidad del lecho m3huecos/m3lecho

# energy
U = 252 #kJ/(m2*h*K)
ts = 300+275.15 #temperatura del aceite termico (K)
t0 = 25+273.15 # temperatura del estado de referencia en K
alpha = np.array([[-1, 1, 0, 0, 1],
                  [ 0,-1, 1, 0, 2],
                  [ 0,-1, 0, 1, 1],
                  [ 0, 0, 1,-1, 1]])

DH0 = np.array([-286.2, -226.1, -96.4, -121.9, 0.0]) #entalpia de cada compoenente a la T de referencia
A = np.array([3.11701732e-10,  4.3670771e-11,  43.400,  87.0774, 9.39393940e-12])
B = np.array([-7.50617201e-7,  6.07692501e-9, 244.500,  258.493, -1.38612795E-8])
C = np.array([ 3.93429762e-4, -3.43061927e-4, 1152.00,  797.321,  6.70727273e-6])
D = np.array([ 3.81447123e-1,  6.31692001e-1, 151.200, -90.5167, -3.83142257e-4])
E = np.array([    1.39426909,     -37.291909, -507.00,  962.718,     28.9049323])

m_aceite = 3. #caudal masico de aceite termico kg/h
Aceite = np.array([2.8645834e-10, -5.05316938e-7, 3.29206037e-4, -9.07074183e-2, 10.909833])

# Constantes
#              [OL, ONA, FENOL, CXENONA, H2]
PM = np.array([100, 98, 94, 96, 2]) # kg/kmol
R = 8.314472/101325*1000 # m3*atm/(K*kmol)
ro_nu = 1.e5 #densidad del gas entre la viscosidad del gas Pa*s
