#    blaance de energia para obtner la variacion de la temperatura
#   con la longitud del reactor
import kinetic2 as kn
import numpy as np
from scipy import integrate


class energia(kn.CuZn):
    '''
    Variacion de a temperatura respecto a la masa decatalizador (dT/dW)
    de cobre zinc los parametros son:
        flujo_molar (kmol/h) = [ciclohexanol, ciclohexanona, fenol, ciclohexenona, H2]
        temperatura (C)
        presion (atm)
        ts (C) temperatura del aceite termico

    Por defecto los valores son:
        flujo_molar (kmol/h) = [0, 0, 0, 0 ,0]
        temperatura (C) = 270
        temperatura del aceite termico (K) = 300
        presion (atm) = 1.0
    '''
    def __init__(self, L, Dint=0.05, Adiabatico=False):

        self.U = 252 #kJ/(m2*h*K)
        self.N = 10 #numero de tubos
        self.L = 10 #longitud de los tubos
        self.ro_l = 2000 #densidad del lecho
        self.L = L
        self.Dint = Dint # diametro interno
        self.ro_l = 1700 #densidad del lecho kg/m3_lecho
        self.ro_p = 2000 #densidad del catalizador kg/m3_cat
        self.e = 1-self.ro_l/self.ro_p #porosidad del lecho m3huecos/m3lecho
        self.Mcat = self.ro_l*self.V_lecho() #masa de catalizados (kg)
        self.ts = 300+275.15 #temperatura del aceite termico (K)
        self.t0 = 25+273.15 # temperatura del estado de referencia en K
        self.alpha = np.array([[-1, 1, 0, 0, 1],
                               [ 0,-1, 1, 0, 2],
                               [ 0,-1, 0, 1, 1],
                               [ 0, 0, 1,-1, 1]])
        self.DH0 = np.array([-286.2, -226.1, -96.4, -121.9, 0.0]) #entalpia de cada compoenente a la T de referencia

             # temperatura en K
             #               Ciclohexanol   Ciclohexanona  Fenol ciclohexenona      H2
        self.A = np.array([3.11701732e-10,  4.3670771e-11,  43.400,  87.0774, 9.39393940e-12])
        self.B = np.array([-7.50617201e-7,  6.07692501e-9, 244.500,  258.493, -1.38612795E-8])
        self.C = np.array([ 3.93429762e-4, -3.43061927e-4, 1152.00,  797.321,  6.70727273e-6])
        self.D = np.array([ 3.81447123e-1,  6.31692001e-1, 151.200, -90.5167, -3.83142257e-4])
        self.E = np.array([    1.39426909,     -37.291909, -507.00,  962.718,     28.9049323])

        self.m_aceite = 500 #caudal masico de aceite termico kg/h
        self.Aceite = np.array([2.8645834e-10, -5.05316938e-7, 3.29206037e-4, -9.07074183e-2, 10.909833])

        self.Adiabatico = Adiabatico

    def V_lecho(self):
        return self.N*self.S(self.Dint)*self.L

    def S(self, D):
        return np.pi*self.Dint**2/4

    def dQdW(self,t, ts):
        if (self.Adiabatico):
            dQdW = 0
        else:
            dQdW = self.U*(np.pi*self.N*self.L/self.ro_l/self.Mcat)**0.5*(ts-t)

        return dQdW

    def Cp(self, t):
        Cp = np.zeros((5,1))

        for i in range(5):
            if i == 0 or i == 1 or  i == 4:
                Cp[i] = self.A[i]*(t)**4+self.B[i]*(t)**3+self.C[i]*(t)**2+self.D[i]*(t)+self.E[i]
            else:
                Cp[i] = self.A[i]+self.B[i]*(self.C[i]/(t)/np.sinh(self.C[i]/(t)))**2+self.D[i]*(self.E[i]/(t)/np.sinh(self.E[i]/(t)))**2

        return Cp

    def Cp1(self, t):
        return self.Cp(t)[0]

    def Cp2(self, t):
        return self.Cp(t)[1]

    def Cp3(self, t):
        return self.Cp(t)[2]

    def Cp4(self, t):
        return self.Cp(t)[3]

    def Cp5(self, t):
        return self.Cp(t)[4]

    def DHi(self, t):

        DHf = np.zeros((5,1)) #entalpia de cada componente
        DHi = np.zeros((4,1)) # entalpia de reaccion

        integracion, _ = integrate.quad(self.Cp1, self.t0, t)
        DHf[0] = self.DH0[0] + integracion
        integracion, _ = integrate.quad(self.Cp2, self.t0, t)
        DHf[1] = self.DH0[1] + integracion
        integracion, _ = integrate.quad(self.Cp3, self.t0, t)
        DHf[2] = self.DH0[2] + integracion
        integracion, _ = integrate.quad(self.Cp4, self.t0, t)
        DHf[3] = self.DH0[3] + integracion
        integracion, _ = integrate.quad(self.Cp5, self.t0, t)
        DHf[4] = self.DH0[4] + integracion

        for i in range(4):
            for j in range(5):
                DHi[i] = self.alpha[i,j]*DHf[j]

        return DHi

    def Cpaceite(self, ts):
        return self.Aceite[0]*(ts)**4+self.Aceite[1]*(ts)**3+self.Aceite[2]*(ts)**2+self.Aceite[3]*(ts)+self.Aceite[4]

    def dTdW(self, n, t, ts, p):
        a = np.dot(self.velocidadReaccion_i(n, t, p), self.DHi(t))
        b = np.dot(n, self.Cp(t))
        return (self.dQdW(t, ts)-a)/b

    def dTsdW(self, n, t, ts, p):
        return self.dTdW(n, t, ts, p)/(self.m_aceite*self.Cpaceite(ts))

    def dTdL(self, n, t, ts, p):
        return self.ro_l*self.N*np.pi*self.Dint**2./4*self.dTdW(n, t, ts, p)

    def dTsdL(self, n, t, ts, p):
        return self.ro_l*self.N*np.pi*self.Dint**2./4*self.dTsdW(n, t, ts, p)




















#
