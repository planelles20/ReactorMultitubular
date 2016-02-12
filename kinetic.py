#   Ecuaciones cineticas de los compuestos quimicos
#   que intervienen en la formacion de ciclohexanona
#   en funcion del catalizador que se use
import numpy as np

class CuZn:
    '''
    Catalizdor de cobre cinc, los parametros a pasar son:
        flujo_molar (kmol/h) = [ciclohexanol, ciclohexanona, fenol, ciclohexenona, H2]
        temperatura (C)
        presion (atm)

        la velocidad para cada componente es en: kmol/(kg_cat*h)

    Por defecto los valores son:
        flujo_molar (kmol/h) = [0, 0, 0, 0 ,0]
        temperatura (C) = 270
        presion (atm) = 1
    '''
    def __init__ (self, Ntub=2400, Dint=0.05):
        self.ro_l = 1700
        self.Dint = Dint
        self.N = Ntub
        self.alpha = np.array([[-1, 1, 0, 0, 1],
                               [ 0,-1, 1, 0, 2],
                               [ 0,-1, 0, 1, 1],
                               [ 0, 0, 1,-1, 1]])
    def y(self, n):
        n = np.array(n)
        return n/sum(n)

    #velocidad de reaccion para el caso A
    def velocidadReaccion_i(self, n, t, p):
        def r1():
            a = np.exp(16.96-6896/t)
            b = (self.y(n)[0]-self.y(n)[4]*self.y(n)[1]/np.exp(13.91-7515/t))*p
            c = (1+0.864*self.y(n)[0]+1.610*self.y(n)[1])*p
            return a*b/c

        def r2():
            a = np.exp(23.65-13506/t)*self.y(n)[1]
            b = (1+0.864*self.y(n)[0]+1.610*self.y(n)[1])*p
            return a/b

        def r3():
            a = np.exp(23.65-13506/t)
            b = (self.y(n)[1]-self.y(n)[4]*self.y(n)[3]/np.exp(11.61-10689/t))*p
            c = (1+0.864*self.y(n)[0]+1.610*self.y(n)[1])*p
            return a*b/c
        def r4():
            a = np.exp(3.51-1298/t)*self.y(n)[3]*p
            b = (1+0.864*self.y(n)[0]+1.610*self.y(n)[1])*p
            return a/b

        return np.array([r1(), r2(), r3(), r4()])/1e3

    def dnjdW(self, n, t, p):
        r_ol = -self.velocidadReaccion_i(n,t,p)[0]
        r_ona = self.velocidadReaccion_i(n,t,p)[0]-self.velocidadReaccion_i(n,t,p)[1]-self.velocidadReaccion_i(n,t,p)[2]
        r_h2 = self.velocidadReaccion_i(n,t,p)[0]+2*self.velocidadReaccion_i(n,t,p)[1]+self.velocidadReaccion_i(n,t,p)[2]+self.velocidadReaccion_i(n,t,p)[3]
        r_cxenona = self.velocidadReaccion_i(n,t,p)[2]-self.velocidadReaccion_i(n,t,p)[3]
        r_phOH = self.velocidadReaccion_i(n,t,p)[1]+self.velocidadReaccion_i(n,t,p)[3]
        return np.array([r_ol, r_ona, r_phOH, r_cxenona, r_h2])

    def dnjdL(self, n, t, p):
        return self.ro_l*self.N*np.pi*self.Dint**2./4*self.dnjdW(n, t, p)









##
