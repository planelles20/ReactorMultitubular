
import numpy as np
import matplotlib.pyplot as plt

class plotearVSCat(object):
    def __init__(self, y, x):
        self.y = y
        self.x = x

    def componentes(self):
        plt.subplot(121)
        plt.plot(self.x,self.y[:,0], label="ciclohexanol")
        plt.plot(self.x,self.y[:,1],label="ciclohexanona")
        plt.plot(self.x,self.y[:,2],label="fenol")
        plt.plot(self.x,self.y[:,3],label="ciclohexenona")
        plt.plot(self.x,self.y[:,4],label="H2")
        plt.xlabel("masa del catalizador (kg)")
        plt.ylabel("flujo molar (kmol/h)")
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.show()

    def T_and_P(self):
        plt.subplot(211)
        plt.plot(self.x,self.y[:,5])
        plt.xlabel("masa del catalizador (kg)")
        plt.ylabel("temperatura (C)")
        plt.subplot(212)
        plt.plot(self.x,self.y[:,7])
        plt.xlabel("masa del catalizador (kg)")
        plt.ylabel("Presion (atm)")
        plt.show()

    def T_and_Ts(self):
        plt.subplot(211)
        plt.plot(self.x,self.y[:,5])
        plt.xlabel("masa del catalizador (kg)")
        plt.ylabel("temperatura (C)")
        plt.subplot(212)
        plt.plot(self.x,self.y[:,6])
        plt.xlabel("masa del catalizador (kg)")
        plt.ylabel("Temperatura del aceite (C)")
        plt.show()

class plotearVSlong(object):
    def __init__(self, y, x):
        self.y = y
        self.x = x

    def componentes(self):
        plt.subplot(121)
        plt.plot(self.x,self.y[:,0], label="ciclohexanol")
        plt.plot(self.x,self.y[:,1],label="ciclohexanona")
        plt.plot(self.x,self.y[:,2],label="fenol")
        plt.plot(self.x,self.y[:,3],label="ciclohexenona")
        plt.plot(self.x,self.y[:,4],label="H2")
        plt.xlabel("longitud del reactor (m)")
        plt.ylabel("flujo molar (kmol/h)")
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.show()

    def T_and_P(self):
        plt.subplot(211)
        plt.plot(self.x,self.y[:,5])
        plt.xlabel("longitud del reactor (m)")
        plt.ylabel("temperatura (C)")
        plt.subplot(212)
        plt.plot(self.x,self.y[:,7])
        plt.xlabel("longitud del reactor (m)")
        plt.ylabel("Presion (atm)")
        plt.show()

    def T_and_Ts(self):
        plt.subplot(211)
        plt.plot(self.x,self.y[:,5])
        plt.xlabel("longitud del reactor (m)")
        plt.ylabel("temperatura (C)")
        plt.subplot(212)
        plt.plot(self.x,self.y[:,6])
        plt.xlabel("longitud del reactor (m)")
        plt.ylabel("Temperatura del aceite (C)")
        plt.show()


class plotearDataVSCat(object):
    def __init__(self, arg):
        self.archivo = str(arg)

    def componentes(self):
        B = np.loadtxt(self.archivo)
        plt.subplot(121)
        plt.plot(B[1:,0],B[1:,1], label="ciclohexanol")
        plt.plot(B[1:,0],B[1:,2],label="ciclohexanona")
        plt.plot(B[1:,0],B[1:,3],label="fenol")
        plt.plot(B[1:,0],B[1:,4],label="ciclohexenona")
        plt.plot(B[1:,0],B[1:,5],label="H2")
        plt.xlabel("Numero de tubos")
        plt.ylabel("flujo molar (kmol/h)")
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.show()

    def T_and_P(self):
        B = np.loadtxt(self.archivo)
        plt.subplot(211)
        plt.plot(B[1:,0],B[1:,6])
        plt.xlabel("Numero de tubos")
        plt.ylabel("temperatura (C)")
        plt.subplot(212)
        plt.plot(B[1:,0],B[1:,8])
        plt.xlabel("Numero de tubos")
        plt.ylabel("Presion (atm)")
        plt.show()

    def T_and_Ts(self):
        B = np.loadtxt(self.archivo)
        plt.subplot(211)
        plt.plot(B[1:,0],B[1:,6])
        plt.xlabel("Numero de tubos")
        plt.ylabel("temperatura (C)")
        plt.subplot(212)
        plt.plot(B[1:,0],B[1:,7])
        plt.xlabel("Numero de tubos")
        plt.ylabel("Temperatura del aceite (C)")
        plt.show()
