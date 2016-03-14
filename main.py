import numpy as np
import matplotlib.pyplot as plt

import representar as rep
import integrateODE as intODE
import estudio as std

if __name__ == "__main__":

    #condiciones iniciales
    n0 = [75.0, 3.75, 0, 0, 0] # kmol/h
    P0 = 2 # atm
    T0 = 270 # C
    Ts0 = 270 # C
    L = 1.5 #m
    N = 2400. #numero de tubos inicial
    Dint = 0.05 # diamtero interno de los tubos (m)

    adi = False

    # que hace?
    estudioNtubosVScat = False
    SolucionParaNTubosVScat = False
    SolucionParaNTubosVSlong = True
    CargarDatosTubos = False

    #Resolver:
    #     n0, T0, Ts, P0      ----->     y = [n0, T0, Ts0, P0]
    # dnj/dW = f(nj, T, P)                       |dnj/dW|
    # dT/dW  = h(nj, T, P)    ----->     dy/dW = | dT/dW|
    #dTs/dW  = i(nj, T, P)                       |dTs/dW|
    # dP/dW  = j(nj, T, P)                       | dP/dW|

    #estudio numero de tubos
    if (estudioNtubosVScat):
        estudio_tubos = std.herraminetas(n0, T0, Ts0, P0)
        estudio_tubos.numeroTubos(N, 0.11, save=adi)

    #solucion para un longitud L y un numero de tubos N
    #catalizador
    if(SolucionParaNTubosVScat): #una solucion
        sol = intODE.ODE(n0, T0, Ts0, P0, L=L, Dint=Dint, Ntub=N, Adiabatico=adi)
        ycat = sol.solutionCat()
        xcat = sol.abcisasMasaCat()
    #longitud
    if(SolucionParaNTubosVSlong): #una solucion
        sol = intODE.ODE(n0, T0, Ts0, P0, L=L, Dint=Dint, Ntub=N, Adiabatico=adi)
        ylong = sol.solutionLong()
        xlong = sol.abcisasLongReactor()


    # representar soluciones
    if(CargarDatosTubos):
        arch = rep.plotearDataVSCat('data/estudioTubos.dat')
        arch.componentes() #
        arch.T_and_P()     #
        arch.T_and_Ts()
    if(SolucionParaNTubosVScat): #una solucion
        rep1 = rep.plotearVSCat(ycat,xcat)
        rep1.componentes() #
        rep1.T_and_P()     #
        rep1.T_and_Ts()
    if(SolucionParaNTubosVSlong): #una solucion
        rep2 = rep.plotearVSlong(ylong,xlong)
        rep2.componentes() #
        rep2.T_and_P()     #
        rep2.T_and_Ts()
