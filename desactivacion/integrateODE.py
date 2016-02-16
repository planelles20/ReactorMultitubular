import numpy as np
from scipy import integrate

import kinetic as kn
import energyBalance as nrgy
import balanceMecanico as mec

class ODE(mec.mecanic):

    def __init__(self, n0, T0, Ts0, P0, a, Dint=0.05, L = 10., Ntub=3000., Adiabatico=False, n = 100):
        # condiciones iniciales
        self.y0 = np.hstack((n0,T0,Ts0,P0))
        self.a = a #es un array
        self.n = n

        #diemsiones de los tubos
        self.L = L  #longitud de los tubos
        self.N = Ntub
        self.Dint = Dint

        self.l = np.linspace(0,L,n)
        self.dl = self.l[1]-self.l[0]

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

        self.www = np.ones((n+1))*5
        self.i = 0


    def funcOdeCat(self, y, x):
        n = y[0:5]
        T = y[5]+273.15
        Ts = y[6]+273.15
        P = y[7]

        # a es un escalar

        dnjdw = self.dnjdW(n, T, P, self.a)
        dTdw = self.dTdW(n, T, Ts, P, self.a)
        dTsdw = self.dTsdW(n, T, Ts, P, self.a)
        dPdw = self.dPdW(n, T, P)
        dydw = np.hstack((dnjdw,dTdw, dTsdw, dPdw))
        return dydw

    def solutionCat(self):
        return  integrate.odeint(self.funcOdeCat, self.y0, self.w)


    def funcOdeLong(self, y, x):
        n = y[0:5]
        T = y[5]+273.15
        Ts = y[6]+273.15
        P = y[7]

        dnjdL = self.dnjdL(n, T, P, self.a)
        dTdL = self.dTdL(n, T, Ts, P, self.a)
        dTsdL = self.dTsdL(n, T, Ts, P, self.a)
        dPdL = self.dPdL(n, T, P)
        dydL = np.hstack((dnjdL,dTdL, dTsdL, dPdL))
        return dydL

    def solutionLong(self):
        return  integrate.odeint(self.funcOdeLong, self.y0, self.l)

    def abcisasMasaCat(self):
        return self.w

    def abcisasLongReactor(self):
        return self.l

    def funcOdeCat2(self, y, x, a):
        n = y[0:5]
        T = y[5]+273.15
        Ts = y[6]+273.15
        P = y[7]

        # a es un escalar

        dnjdw = self.dnjdW(n, T, P, a)
        dTdw = self.dTdW(n, T, Ts, P, a)
        dTsdw = self.dTsdW(n, T, Ts, P, a)
        dPdw = self.dPdW(n, T, P)
        dydw = np.hstack((dnjdw,dTdw, dTsdw, dPdw))
        return dydw


    def funcOdeLong2(self, y, x, a):
        n = y[0:5]
        T = y[5]+273.15
        Ts = y[6]+273.15
        P = y[7]

        dnjdL = self.dnjdL(n, T, P, a)
        dTdL = self.dTdL(n, T, Ts, P, a)
        dTsdL = self.dTsdL(n, T, Ts, P, a)
        dPdL = self.dPdL(n, T, P)
        dydL = np.hstack((dnjdL,dTdL, dTsdL, dPdL))
        return dydL

    def Runge_Kutta(self, f, y0, x):

        ysol = np.zeros((x.size, y0.size))
        ysol[0,:] = y0
        h = x[1]-x[0]

        for i in range(1,self.n):
            a = self.a[i]
            k1 = f(ysol[i-1,:]        , x     , a)
            k2 = f(ysol[i-1,:]+k1*h/2., x+h/2., a)
            k3 = f(ysol[i-1,:]+k2*h/2., x+h/2., a)
            k4 = f(ysol[i-1,:]+k3*h   , x+h   , a)

            ysol[i,:] = ysol[i-1,:]+h/6.*(k1+2*k2+2*k3+k4)

        return ysol

    def solutionLong2(self):
        return  self.Runge_Kutta(self.funcOdeLong2, self.y0, self.l)


    def aaa(self, nj0, T, P, a, dt):
        '''
            devuelve un array de variacion de la actividad
            nj0 es un array de (nx5)
            T es un array de (nx1)
            P es un array de (nx1)
            a es un array de (nx1)
        '''
        da = np.zeros((a.size), dtype=float)
        for i in range(a.size):
            a1 = 0.036*self.y(nj0[i,:])[0]*P[i]+0.069*self.y(nj0[i,:])[1]*P[i]
            a2 = 1+2.132*self.y(nj0[i,:])[4]*P[i]
            a3 = a[i]-0.408/(1+1.975*self.y(nj0[i,:])[1]*P[i])
            da[i] = -a1/a2*a3

        return a+da*dt







        #
