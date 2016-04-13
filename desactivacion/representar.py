
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import datos
import kinetic as kn


class plotearLong(object):
    def __init__(self, y, x, t):
        self.SOL = y
        xx, tt = np.meshgrid(x, t)
        self.X = xx
        self.Y = tt

        self.nl = len(x)
        self.nt = len(t)

    def componentes(self):

        fig = plt.figure(figsize=plt.figaspect(0.5))

        ax = fig.add_subplot(2, 3, 1, projection='3d')
        surf = ax.plot_surface(self.X, self.Y,self.SOL[:,:,0],cmap=cm.coolwarm)
        plt.xlabel("longitud (m)")
        plt.ylabel("tiempo (horas)")
        plt.title("ciclohexanol (kmol/h)")

        ax = fig.add_subplot(2, 3, 2, projection='3d')
        surf = ax.plot_surface(self.X, self.Y, self.SOL[:,:,1],cmap=cm.coolwarm)
        plt.xlabel("longitud (m)")
        plt.ylabel("tiempo (horas)")
        plt.title("ciclohexanona (kmol/h)")

        ax = fig.add_subplot(2, 3, 3, projection='3d')
        surf = ax.plot_surface(self.X, self.Y, self.SOL[:,:,2],cmap=cm.coolwarm)
        plt.xlabel("longitud (m)")
        plt.ylabel("tiempo (horas)")
        plt.title("fenol (kmol/h)")

        ax = fig.add_subplot(2, 3, 4, projection='3d')
        surf = ax.plot_surface(self.X, self.Y, self.SOL[:,:,3]*1000,cmap=cm.coolwarm)
        plt.xlabel("longitud (m)")
        plt.ylabel("tiempo (horas)")
        plt.title("ciclohexenona (mol/h)")

        ax = fig.add_subplot(2, 3, 6, projection='3d')
        surf = ax.plot_surface(self.X, self.Y, self.SOL[:,:,4],cmap=cm.coolwarm)
        plt.xlabel("longitud (m)")
        plt.ylabel("tiempo (horas)")
        plt.title("H2 (kmol/h)")

        plt.show()


    def T_and_P(self):
        fig = plt.figure(figsize=plt.figaspect(0.5))

        ax = fig.add_subplot(1, 2, 1, projection='3d')
        surf = ax.plot_surface(self.X, self.Y,self.SOL[:,:,5],cmap=cm.coolwarm)
        plt.xlabel("longitud (m)")
        plt.ylabel("tiempo (horas)")
        plt.title("Temperatura del fluido(C)")

        ax = fig.add_subplot(1, 2, 2, projection='3d')
        surf = ax.plot_surface(self.X, self.Y,self.SOL[:,:,7],cmap=cm.coolwarm)
        plt.xlabel("longitud (m)")
        plt.ylabel("tiempo (horas)")
        plt.title("Presion(atm)")

        plt.show()

    def T_and_Ts(self):
        fig = plt.figure(figsize=plt.figaspect(0.5))

        ax = fig.add_subplot(1, 2, 1, projection='3d')
        surf = ax.plot_surface(self.X, self.Y,self.SOL[:,:,5],cmap=cm.coolwarm)
        plt.xlabel("longitud (m)")
        plt.ylabel("tiempo (horas)")
        plt.title("Temperatura del fluido(C)")

        ax = fig.add_subplot(1, 2, 2, projection='3d')
        surf = ax.plot_surface(self.X, self.Y,self.SOL[:,:,6],cmap=cm.coolwarm)
        plt.xlabel("longitud (m)")
        plt.ylabel("tiempo (horas)")
        plt.title("Temperatura del refrigerante (C)")

        plt.show()

    def actividad(self):
        fig = plt.figure(figsize=plt.figaspect(0.5))

        ax = fig.add_subplot(1, 1, 1, projection='3d')
        surf = ax.plot_surface(self.X, self.Y,self.SOL[:,:,8],cmap=cm.coolwarm)
        plt.xlabel("longitud (m)")
        plt.ylabel("tiempo (horas)")
        plt.title("actividad")

        plt.show()

    def moduloThiele(self):
        mL = np.zeros((self.nt,self.nl))
        rr = kn.CuZn()
        for i in range(self.nt):
            for j in range(self.nl):
                n = self.SOL[i,j,0:5] #kmol/h
                T = self.SOL[i,j,5]+273.15# K
                P = self.SOL[i,j,7] #atm
                a = self.SOL[i,j,8]
                L = datos.Dint/2/2 #m
                r =  -rr.dnjdW(n,T,P,a)[0]/3600/1000 #kmol/(kgs), (tiene que ser positiva)
                if (r<0): r = 0
                C_ol = n[0]/((sum(n))*datos.R*T/P) #kmol/m3
                k = r*datos.ro_l/C_ol # s(-1)
                D_OL = datos.D0_OL*np.exp(-datos.Ed_OL/((datos.R*101325)*T)) #cambiar la R a kJ/molK
                De = D_OL*datos.e_cat**2 #m2/s
                X = 0.712505333333 #grado de conversion de ciclohexanol en equilibrio
                mL[i,j] = L*(k/De/X)**0.5 #adimensional

        fig = plt.figure(figsize=plt.figaspect(0.5))
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        surf = ax.plot_surface(self.X, self.Y, mL, cmap=cm.coolwarm)
        plt.xlabel("longitud (m)")
        plt.ylabel("tiempo (horas)")
        plt.title("modulo de Thiele")
        plt.show()


    def general_plot(self, y, x, t, title=" "):
        X, Y = np.meshgrid(x, t)
        fig = plt.figure(figsize=plt.figaspect(0.5))
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        surf = ax.plot_surface(X, Y, y,cmap=cm.coolwarm)
        plt.xlabel("longitud (m)")
        plt.ylabel("tiempo (horas)")
        plt.title(title)
        plt.show()
