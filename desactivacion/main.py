import numpy as np
import matplotlib.pyplot as plt

import representar as rep
import integrateODE as intODE

if __name__ == "__main__":

    #condiciones iniciales
    n0 = [75.0, 3.75, 0, 0, 0] # kmol/h
    P0 = 2 # atm
    T0 = 270 # C
    Ts0 = 270 # C
    L = 2.5 #m
    N = 2400. #numero de tubos inicial
    Dint = 0.05 # diamtero interno de los tubos (m)

    adi = False

    nl = 100
    nt = 20
    tf = 150 # horas

    # que hace?

    #Resolver:
    #     n0, T0, Ts, P0      ----->     y = [n0, T0, Ts0, P0]
    # dnj/dW = f(nj, T, P, a)                    |dnj/dW|
    # dT/dW  = h(nj, T, P, a) ----->     dy/dW = | dT/dW|
    #dTs/dW  = i(nj, T, P, a)                    |dTs/dW|
    # dP/dW  = j(nj, T, P)                       | dP/dW|
    # dadt   = k(nj, T, P, a)                    | da/dt|

    #solucion para un longitud L y un numero de tubos N
    #catalizador
    if(False): #una solucion
        sol = intODE.ODE(n0, T0, Ts0, P0, a, L=L, Dint=Dint, Ntub=N, Adiabatico=adi)
        ycat = sol.solutionCat()
        xcat = sol.abcisasMasaCat()

    #longitud
    SOL = np.zeros((nt,nl,9))
    tlong = np.linspace(0,tf,nt)
    dt = tlong[1]-tlong[0]
    a = np.ones((nl))
    for i in range(nt):
        sol = intODE.ODE(n0, T0, Ts0, P0, a, L=L, Dint=Dint, Ntub=N, Adiabatico=adi, n=nl)
        ylong = sol.solutionLong2()
        SOL[i,:,:8] = ylong
        SOL[i,:,8] = a
        #                nj0       T            P            a
        a = sol.aaa(SOL[i,:,:5], SOL[i,:,5], SOL[i,:,7], SOL[i,:,8], dt)

    xlong = sol.abcisasLongReactor()


    rep1 = rep.plotearLong(SOL,xlong,tlong)
    rep1.componentes()
    rep1.T_and_P()     #
    rep1.T_and_Ts()
    rep1.actividad()
    # representar soluciones

    plt.plot(tlong, SOL[:,-1,:5])
    plt.show()

    if(False): #una solucion
        rep1 = rep.plotearVSCat(ycat,xcat)
        rep1.componentes() #
        rep1.T_and_P()     #
        rep1.T_and_Ts()
