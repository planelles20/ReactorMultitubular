import numpy as np
import datos
import matplotlib.pyplot as plt

#datos para la dispersion
u = 5 # velocidad del viento (m/s)
te = 10*60 #tiempo de emision de la fuga (segundos)
y = 0 #ordenada y (metros)
z = 0 #ordenada z (metrosS)
h = 0 #altura de la fuga

#coeficientes de dispersion
a = 0.128
b = 0.905
c = 0.20
d = 0.76

###########
nl = datos.nl
nt = datos.nt
X0 = 0.05
X = 1.5
nx = 100
x = np.linspace(X0,X,nx) #distancia del eje x en metros

# limites de inflamabilidad en % de volumen
LIIe = 1*np.ones((nx)) #limite infereior de inflamabilidad de la mezcla a la entrada al reactor
LIIs = 1.5*np.ones((nx))#limite infereior de inflamabilidad de la mezcla a la salida al reactor
LSIe = 9*np.ones((nx))#limite superior de inflamabilidad de la mezcla a la entrada al reactor
LSIs = 13.9*np.ones((nx))#limite superior de inflamabilidad de la mezcla a la salida al reactor

N = np.zeros((nt,nl,5))
T = np.zeros((nt,nl,1))
P = np.zeros((nt,nl,1))
Ts = np.zeros((nt,nl,1))

N[:,:,0] = np.loadtxt('./data/desactivacion0.dat')
N[:,:,1] = np.loadtxt('./data/desactivacion1.dat')
N[:,:,2] = np.loadtxt('./data/desactivacion2.dat')
N[:,:,3] = np.loadtxt('./data/desactivacion3.dat')
N[:,:,4] = np.loadtxt('./data/desactivacion4.dat')
T[:,:,0] = np.loadtxt('./data/desactivacion5.dat')
P[:,:,0] = np.loadtxt('./data/desactivacion6.dat')

Ne = N[0,0,:]  #flujo molar a la entrada del reactor (kmol/h)
Ns = N[0,-1,:] #flujo molar a la salida del reactor (kmol/h)
Te = T[0,0,0]  #temperatura a la entrada del reactor (C)
Ts = T[0,-1,0] #temperatura a la salida del reactor (C)
Pe = P[0,0,0]  #presion a la entrada del reactor (atm)
Ps = P[0,-1,0] #presion a la salida del reactor (atm)

#distncia entre instantanea y continua
x_ic = 1.8*u*te

# calculo de la masa de emision a la entrada del reactor (Me) y a la salida (Ms)
Qe = Ne*datos.PM/3600
Qs = Ns*datos.PM/3600
Q = sum(Ns*datos.PM)/3600
# calculo del perfil de concentracion en el eje x
Conc_e = np.zeros((nx,5)) #kg/m3 entrada al reactor por componentes
Conc_s = np.zeros((nx,5)) #kg/m3 salida al reactor por componentes
Conc = np.zeros((nx)) #kg/m3 salida al reactor en total
V_e = np.zeros((nx)) # % volumen de aire de la mezcla a la entrada al reactor
V_s =  np.zeros((nx)) # % volumen de aire de la mezcla a la salida al reactor

#volumen inical que ocupan los reactivos al salir por la fuga
V_e0 = None
V_s0 = None

for i in range(nx):
    if x[i]<x_ic:
        sx = 0.13*x[i]
        sy = 0.5*a*x[i]**b
        sz = c*x[i]**d
    else:
        sx = x[i]
        sy = a*x[i]**b
        sz = c*x[i]**d

    for j in range(5):
        Conc_e[i,j] = Qe[j]/(np.pi*sy*sz*u)*np.exp(-0.5*(y**2/sy**2+z**2/sz**2))
        Conc_s[i,j] = Qs[j]/(np.pi*sy*sz*u)*np.exp(-0.5*(y**2/sy**2+z**2/sz**2))

    Conc[i] = Q/(np.pi*sy*sz*u)*np.exp(-0.5*(y**2/sy**2+z**2/sz**2))

    #calculo en % en volumen
    V_e[i] = sum(Ne)*datos.R*(Te+273.15)/Pe/(np.pi*sy*sz*u)*np.exp(-0.5*(y**2/sy**2+z**2/sz**2))
    V_s[i] = sum(Ns)*datos.R*(Ts+273.15)/Ps/(np.pi*sy*sz*u)*np.exp(-0.5*(y**2/sy**2+z**2/sz**2))

    if i == 0:
        V_e0 = V_e[i]
        V_s0 = V_s[i]

    V_e[i] = V_e[i]/V_e0*100
    V_s[i] = V_s[i]/V_s0*100


#Representar
fig1 = plt.figure(1)
plt.title("Dispersion de los componentes a la entrada")
plt.xlabel("longitud (m)")
plt.xlim([X0,X])
plt.ylabel("Concentracion (kg/m3)")
for j in range(5):
    plt.plot(x, Conc_e[:,j])
plt.show()

fig3 = plt.figure(1)
plt.title("Dispersion de los componentes a la salida")
plt.xlabel("longitud (m)")
plt.xlim([X0,X])
plt.ylabel("Concentracion (kg/m3)")
for j in range(5):
    plt.plot(x, Conc_s[:,j])
plt.show()


fig2 = plt.figure(2)
plt.title("Dispersion caudal total")
plt.xlabel("longitud (m)")
plt.xlim([X0,X])
plt.ylabel("Concentracion (kg/m3)")
plt.plot(x, Conc)
plt.show()

fig4 = plt.figure(2)
plt.title("Dispersion caudal total a la entrada")
plt.xlabel("longitud (m)")
plt.xlim([X0,X])
plt.ylabel("Porcentaje en volumen")
plt.plot(x, V_e)
plt.plot(x, LIIe)
plt.plot(x, LSIe)
plt.show()

fig4 = plt.figure(2)
plt.title("Dispersion caudal total a la salida")
plt.xlabel("longitud (m)")
plt.xlim([X0,X])
plt.ylabel("Porcentaje en volumen")
plt.plot(x, V_s)
plt.plot(x, LIIs)
plt.plot(x, LSIs)
plt.show()
