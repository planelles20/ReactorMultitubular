import numpy as np

# diseño del intercambiador
Dint = 2.*2.54/100 # metros, diametro interno de los tubos
Ntub = 3000 # nuemro de tubos
L = 2.5 # metros, longitud de los tubos



# catalizador
ro_l = 1760. #densidad del lecho kg/m3_lecho
ro_p = 2300. #densidad del catalizador kg/m3_cat
dp = 3e-3 #diametro medio del Catalizdor (m)
e = 1-ro_l/ro_p #porosidad del lecho m3huecos/m3lecho

# datos para difusion en el catalizador
D0_OL = 7.16e-6 # difusividad del ciclohexanol m2/s (preexponencial factor)
D0_ONA = 7.55e-7 # difusividad del ciclohexanona m2/s (preexponencial factor)
Ed_OL = 31.98 # energia de activacion de la difusion para el ciclohexanol (kJ/mol)
Ed_ONA = 24.25 # energia de activacion de la difusion para el ciclohexanona (kJ/mol)
e_cat = 0.60 # huecos libres dentro del catalizador

# energy
U = 12*5.6745/1000*3600 #kJ/(m2*h*K)
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

#              [OL, ONA, FENOL, CXENONA, H2]
PM = np.array([100, 98, 94, 96, 2]) # kg/kmol
R = 8.314472/101325*1000 # m3*atm/(K*kmol)
ro_nu = 1.e5 #densidad del gas entre la viscosidad del gas kg/(m3*Pa*s)
viscosidad = 1e-5 #viscosidad del gas en Pa s
