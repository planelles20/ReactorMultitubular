import numpy as np
import integrateODE as intODE


class herraminetas(object):
    """docstring for """
    def __init__(self, n0, T0, Ts0, P0, L=10000):
        self.n0 = n0
        self.T0 = T0
        self.Ts0 = Ts0
        self.P0 = P0
        self.L = L


    def numeroTubos(self, N0, caida_presion, dN=10, save=False):
        '''
            (caida_presion > 0.11)
        '''
        N = N0
        y = np.ones((1,8))*2
        A = np.zeros((1,9))
        while y[-1, 7]>caida_presion*self.P0:
            sol = intODE.ODE(self.n0, self.T0, self.Ts0, self.P0, L=self.L, Ntub=N)
            y = sol.solutionCat()
            b = np.hstack((N,y[-1,0],y[-1,1],y[-1,2],y[-1,3],y[-1,4],y[-1,5],y[-1,6],y[-1,7]))
            A = np.vstack((A,b))
            N -= dN

        if (save):
            np.savetxt('data/estudioTubosCat.dat', A, fmt='%.5e')

        return 0

    def longitudReactor(self, N, L):
        #TODO
        pass

    def tubosANDDint(self, N, Dint):
        #TODO
        pass
