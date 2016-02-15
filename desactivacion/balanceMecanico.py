import numpy as np

import energyBalance as nrgy


class mecanic(nrgy.energia):

    def __init__(self, L=10., Ntub=3000, Dint = 0.05):

        #              [OL, ONA, FENOL, CXENONA, H2]
        self.PM = np.array([100, 98, 94, 96, 2]) # kg/kmol
        self.R = 8.3144/101325 # kmol/(K*atm)
        self.N = 2000 #numero de tubos
        self.Dint =  Dint #diametro interno de los tubos (m)
        self.L = L #longitud de los tubos (m)
        self.dp = 3e-3 #diametro medio del Catalizdor (m)
        self.ro_nu = 1e5/3600 #densidad del gas entre la viscosidad del gas Pa*h
        self.ro_l = 1700 #densidad del lecho kg/m3_lecho
        self.ro_p = 2000 #densidad del catalizador kg/m3_cat
        self.e = 1-self.ro_l/self.ro_p #porosidad del lecho m3huecos/m3lecho
        self.Mcat = self.ro_l*self.V_lecho()
        self.N = Ntub

    def V_lecho(self):
        return self.N*self.S(self.Dint)*self.L

    def PMm(self, n):
        return np.dot(self.y(n), self.PM)

    def ro(self, n, t, p):
        #ley gases ideales
        return p*self.PMm(n)/(self.R*t)

    def Qv(self, n, t, p):
        #ley gases ideales
        return sum(n)*self.R*t/p

    def S(self, D):
        return np.pi*D**2/4

    def velocidadGas(self, n, t, p):
        return self.Qv(n, t, p)/(self.N*self.S(self.D))

    def dPdW(self, n, t, p):
        viscosidad = (1/self.ro_nu)*self.ro(n,t,p)
        a = 150*(1-self.e)**2*viscosidad*self.Qv(n,t,p)/3600
        b = self.ro_l*(self.S(self.Dint)*self.N)**2*self.e**3*self.dp**2

        c = 1.75*(1-self.e)*self.ro(n,t,p)*(self.Qv(n,t,p)/3600)**2
        d = self.ro_l*(self.S(self.Dint)*self.N)**3*self.e**3*self.dp
        return -(a/b + c/d)/101325

    def dPdL(self, n, t, p):
        return self.ro_l*self.N*np.pi*self.Dint**2/4*self.dPdW(n, t, p)
