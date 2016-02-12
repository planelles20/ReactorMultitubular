import numpy as np
from scipy import integrate

import kinetic as kn
import energyBalance as nrgy
import balanceMecanico as mec

class ODE(mec.mecanic):

    def __init__(self, n0, T0, Ts0, P0, Dint=0.05, L = 10., Ntub=3000., Adiabatico=False, n = 10000):
        # condiciones iniciales
        self.y0 = np.hstack((n0,T0,Ts0,P0))

        #diemsiones de los tubos
        self.L = L  #longitud de los tubos
        self.N = Ntub
        self.Dint = Dint

        self.l = np.linspace(0,L,n)

        # catalizador
        self.ro_l = 1700. #densidad del lecho kg/m3_lecho
        self.ro_p = 2000. #densidad del catalizador kg/m3_cat
        self.dp = 3.e-3 #diametro medio del Catalizdor (m)
        self.e = 1-self.ro_l/self.ro_p #porosidad del lecho m3huecos/m3lecho
        self.Mcat = self.ro_l*self.V_lecho()

        self.w = np.linspace(0, self.Mcat, n)

        # energy
        self.U = 252 #kJ/(m2*h*K)
        self.ts = 300+275.15 #temperatura del aceite termico (K)
        self.t0 = 25+273.15 # temperatura del estado de referencia en K
        self.alpha = np.array([[-1, 1, 0, 0, 1],
                               [ 0,-1, 1, 0, 2],
                               [ 0,-1, 0, 1, 1],
                               [ 0, 0, 1,-1, 1]])

        self.DH0 = np.array([-286.2, -226.1, -96.4, -121.9, 0.0]) #entalpia de cada compoenente a la T de referencia
        self.A = np.array([3.11701732e-10,  4.3670771e-11,  43.400,  87.0774, 9.39393940e-12])
        self.B = np.array([-7.50617201e-7,  6.07692501e-9, 244.500,  258.493, -1.38612795E-8])
        self.C = np.array([ 3.93429762e-4, -3.43061927e-4, 1152.00,  797.321,  6.70727273e-6])
        self.D = np.array([ 3.81447123e-1,  6.31692001e-1, 151.200, -90.5167, -3.83142257e-4])
        self.E = np.array([    1.39426909,     -37.291909, -507.00,  962.718,     28.9049323])

        self.m_aceite = 5000 #caudal masico de aceite termico kg/h
        self.Aceite = np.array([2.8645834e-10, -5.05316938e-7, 3.29206037e-4, -9.07074183e-2, 10.909833])

        # Constantes
        #              [OL, ONA, FENOL, CXENONA, H2]
        self.PM = np.array([100, 98, 94, 96, 2]) # kg/kmol
        self.R = 8.314472/101325*1000 # m3*atm/(K*kmol)
        self.ro_nu = 1.e5 #densidad del gas entre la viscosidad del gas Pa*s

        self.Adiabatico = Adiabatico

    def funcOdeCat(self, y, w):
        n = y[0:5]
        T = y[5]+273.15
        Ts = y[6]+273.15
        P = y[7]

        dnjdw = self.dnjdW(n, T, P)
        dTdw = self.dTdW(n, T, Ts, P)
        dTsdw = self.dTsdW(n, T, Ts, P)
        dPdw = self.dPdW(n, T, P)
        dydw = np.hstack((dnjdw,dTdw, dTsdw, dPdw))
        return dydw

    def solutionCat(self):
        return  integrate.odeint(self.funcOdeCat, self.y0, self.w)


    def funcOdeLong(self, y, w):
        n = y[0:5]
        T = y[5]+273.15
        Ts = y[6]+273.15
        P = y[7]

        dnjdL = self.dnjdL(n, T, P)
        dTdL = self.dTdL(n, T, Ts, P)
        dTsdL = self.dTsdL(n, T, Ts, P)
        dPdL = self.dPdL(n, T, P)
        dydL = np.hstack((dnjdL,dTdL, dTsdL, dPdL))
        return dydL

    def solutionLong(self):
        return  integrate.odeint(self.funcOdeLong, self.y0, self.l)



    def abcisasMasaCat(self):
        return self.w

    def abcisasLongReactor(self):
        return self.l
