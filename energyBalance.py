#    blaance de energia para obtner la variacion de la temperatura
#   con la longitud del reactor
import kinetic as kn
import numpy as np
from scipy import integrate

import datos

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
    def __init__(self, Adiabatico=False):

        self.U = datos.U #kJ/(m2*h*K)
        self.N = datos.Ntub #numero de tubos
        self.L = datos.L #longitud de los tubos
        self.ro_l = datos.ro_l #densidad del lecho
        self.L = datos.L
        self.Dint = datos.Dint # diametro interno
        self.ro_l = datos.ro_l #densidad del lecho kg/m3_lecho
        self.ro_p = ratos.ro_p #densidad del catalizador kg/m3_cat
        self.e =  datos.e #porosidad del lecho m3huecos/m3lecho
        self.Mcat = self.ro_l*self.V_lecho() #masa de catalizados (kg)
        self.ts = datos.ts #temperatura del aceite termico (K)
        self.t0 = datos.t0 # temperatura del estado de referencia en K
        self.alpha = datos.alpha
        self.DH0 = datos.DH0 #entalpia de cada compoenente a la T de referencia

             # temperatura en K
        self.A = datos.A
        self.B = datos.B
        self.C = datos.C
        self.D = datos.D
        self.E = datos.E

        self.Adiabatico = Adiabatico

    def V_lecho(self):
        return self.N*self.S(self.Dint)*self.L

    def S(self, D):
        return np.pi*self.Dint**2/4

    def dQdW(self,t, ts):
        if (self.Adiabatico):
            dQdW = 0
        else:
            dQdW = 4*self.U*(ts-t)/(self.ro_p*(1-self.e)*self.Dint)
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

    def dTdW(self, n, t, ts, p, a):
        aa = np.dot(self.velocidadReaccion_i(n, t, p, a), self.DHi(t))
        b = np.dot(n, self.Cp(t))
        return (self.dQdW(t,ts)-aa)/b

    def dTsdW(self, n, t, ts, p, a):
        return 0

    def dTdL(self, n, t, ts, p, a):
        return self.ro_l*self.N*np.pi*self.Dint**2./4*self.dTdW(n, t, ts, p, a)

    def dTsdL(self, n, t, ts, p, a):
        return 0



















#
