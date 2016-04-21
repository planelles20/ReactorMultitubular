import numpy as np
import matplotlib.pyplot as plt

import representar as rep
import integrateODE as intODE
import datos
import kinetic as kn

if __name__ == "__main__":

    #condiciones iniciales
    n0 = [75.0, 3.75, 0, 0, 0] # kmol/h
    P0 = 2 # atm
    T0 = 270 # C
    Ts0 = 270 # C

    L = datos.L # metros
    N = datos.Ntub #numero de tubos inicial
    Dint = datos.Dint # diamtero interno de los tubos (m)

    adi = False

    nl = 100 # puntos en la lungitud L
    nt = 8 #puntos en el tiempo
    tf = 50 # horas

    # que hace?

    #Resolver:
    #     n0, T0, Ts, P0      ----->     y = [n0, T0, Ts0, P0, a0]
    # dnj/dW = f(nj, T, P, a)                    |dnj/dW|
    # dT/dW  = h(nj, T, P, a) ----->     dy/dW = | dT/dW|
    #dTs/dW  = i(nj, T, P, a)                    |dTs/dW|
    # dP/dW  = j(nj, T, P)                       | dP/dW|
    # dadt   = k(nj, T, P, a)                    | da/dt|

    #solucion para un longitud L y un numero de tubos N
    #catalizador
    if(False): #una solucion
        sol = intODE.ODE(n0, T0, Ts0, P0, 1, Adiabatico=adi)
        ycat = sol.solutionCat()
        xcat = sol.abcisasMasaCat()

        rep0 = rep.plotearLong(sol,xcat,tlong)
        rep0.componentes()

    #longitud
    if(False):
        SOL = np.zeros((nt,nl,9))
        tlong = np.linspace(0,tf,nt)
        dt = tlong[1]-tlong[0]
        a = np.ones((nl))
        for i in range(nt):
            print (nt-i)
            sol = intODE.ODE(n0, T0, Ts0, P0, a, Adiabatico=adi, n=nl)
            ylong = sol.solutionLong2()
            SOL[i,:,:8] = ylong
            SOL[i,:,8] = a
            #                nj0       T            P            a
            a = sol.aaa(SOL[i,:,:5], SOL[i,:,5], SOL[i,:,7], SOL[i,:,8], dt)

        xlong = sol.abcisasLongReactor()
        #print (SOL[-1,-1,:])
        np.savetxt('../data/desactivacionX.dat', xlong, fmt='%.5e')
        np.savetxt('../data/desactivacion0.dat', SOL[:,:,0], fmt='%.5e')
        np.savetxt('../data/desactivacion1.dat', SOL[:,:,1], fmt='%.5e')
        np.savetxt('../data/desactivacion2.dat', SOL[:,:,2], fmt='%.5e')
        np.savetxt('../data/desactivacion3.dat', SOL[:,:,3], fmt='%.5e')
        np.savetxt('../data/desactivacion4.dat', SOL[:,:,4], fmt='%.5e')
        np.savetxt('../data/desactivacion5.dat', SOL[:,:,5], fmt='%.5e')
        np.savetxt('../data/desactivacion6.dat', SOL[:,:,6], fmt='%.15e')
        np.savetxt('../data/desactivacion7.dat', SOL[:,:,7], fmt='%.5e')
        np.savetxt('../data/desactivacion8.dat', SOL[:,:,8], fmt='%.5e')

    if (True):
        tlong = np.linspace(0,tf,nt)
        xlong = np.zeros((nl))
        xlong = np.loadtxt('../data/desactivacionX.dat')
        SOL = np.zeros((nt,nl,9))
        SOL[:,:,0] = np.loadtxt('../data/desactivacion0.dat')
        SOL[:,:,1] = np.loadtxt('../data/desactivacion1.dat')
        SOL[:,:,2] = np.loadtxt('../data/desactivacion2.dat')
        SOL[:,:,3] = np.loadtxt('../data/desactivacion3.dat')
        SOL[:,:,4] = np.loadtxt('../data/desactivacion4.dat')
        SOL[:,:,5] = np.loadtxt('../data/desactivacion5.dat')
        SOL[:,:,6] = np.loadtxt('../data/desactivacion6.dat')
        SOL[:,:,7] = np.loadtxt('../data/desactivacion7.dat')
        SOL[:,:,8] = np.loadtxt('../data/desactivacion8.dat')
        '''
        rep1 = rep.plotearLong(SOL,xlong,tlong)
        rep1.componentes()
        rep1.T_and_P()
        rep1.T_and_Ts()
        rep1.actividad()
        rep1.general_plot(SOL[:,:,0], xlong, tlong, title="ciclohexanol (kmol/h)")
        rep1.general_plot(SOL[:,:,1], xlong, tlong, title="ciclohexanona (kmol/h)")
        rep1.general_plot(SOL[:,:,5], xlong, tlong, title="Temperatura (C)")
        rep1.general_plot(SOL[:,:,7], xlong, tlong, title="Presion (atm)")
        rep1.moduloThiele()
        rep1.fraccionMolar(n=0,title="ciclohexanol")
        rep1.fraccionMolar(n=1,title="ciclohexanona")
        rep1.conversion
        '''
        X1 = np.zeros((nl))
        X2 = np.zeros((nl))
        X3 = np.zeros((nl))
        X4 = np.zeros((nl))
        X5 = np.zeros((nl))

        for i in range(nl):
            X1[i] = SOL[0,i,0]/sum(SOL[0,i,0:5])
            X2[i] = SOL[0,i,1]/sum(SOL[0,i,0:5])
            X3[i] = SOL[0,i,2]/sum(SOL[0,i,0:5])
            X4[i] = SOL[0,i,3]/sum(SOL[0,i,0:5])
            X5[i] = SOL[0,i,4]/sum(SOL[0,i,0:5])
        # representar
        fig = plt.figure(1)
        plt.title("Flaccion molar")
        plt.xlabel("longitud (m)")
        plt.ylabel("Fraccion molar")
        plt.plot(xlong, X1)
        plt.plot(xlong, X2)
        plt.plot(xlong, X3)
        plt.plot(xlong, X4)
        plt.plot(xlong, X5)
        plt.show()
