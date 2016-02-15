
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


class plotearLong(object):
    def __init__(self, y, x, t):
        self.SOL = y
        xx, tt = np.meshgrid(x, t)
        self.X = xx
        self.Y = tt

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
        surf = ax.plot_surface(self.X, self.Y, self.SOL[:,:,3],cmap=cm.coolwarm)
        plt.xlabel("longitud (m)")
        plt.ylabel("tiempo (horas)")
        plt.title("ciclohexenona (kmol/h)")

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
