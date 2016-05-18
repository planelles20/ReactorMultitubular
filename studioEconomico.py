import numpy as np
import datos

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

#datos
D = 3.64 #metros
Sreactor = np.pi*(D/2)**2 #seccion del reactor m2
Ntub = datos.Ntub
Stub = np.pi*(datos.Dint/2)**2 #seccion de los tubos m2
Cir_tub = np.pi*datos.Dint  #circunferencia de los tubos m

#indices de Marshall y Swif
MS2015 = 1600
MS2001 = 1094

#cambio de dolares a euros
f = 0.876

#conversio de metros a pies
m_ft = 1/3.28084 #m/ft

#tipo de material empleado
coef1 = 1.8 #acero inoxidable
coef2 = 1.0 #acero

# precios de los componentes
PriceONA = 1800/1000 # euros/kg
PriceOL = 1150/1000 # euros/kg
PriceAceite = 2364.52 # euros/m3
PriceCat = 200 #euros/kg

#tiempo de amortizacion
amort = 10 #anos

#horas trabajadas a la semana
horas_semana = 168

#peso molecular
PM_OL = datos.PM[0]   # kmol/kg
PM_ONA = datos.PM[1]  # kmol/kg

#Presion de diseno
P = 1*(101325/1e5)*1.10 #bares de presion

# Vlores de las correlaciones
K1 = 3.5565
K2 = 0.3776
K3 = 0.0905
Fp = ((P+1)*D/(2*(850-0.6*(P+1))+0.00315))/0.0063 #factor de la presion
B1 = 1.49
B2 = 1.52
Fm = 2.8 #factor material

#ecuacion calcula el coste en funcion de la composicion, longitud, tiempo
# unidades: euros/ano
def costeReactor(n_OL, n_ONA, L, t):
    Vtub = Stub*Ntub*L #volumen de los tubos
    Vreac = Sreactor*L-Vtub #volumen que acua el aceite
    A = Cir_tub*Ntub*L #area intercambio de calor
    V = Sreactor*L+1e-10
    #horas que rabaja el rector al ano, teniendo en cuenta que necesita 8 horas
    #para regenerar  el catalizador
    if t > 8:
        Horas_Ano = 12*4*(horas_semana/t)*(t-8)
    else:
        Horas_Ano = 0

    cp = 10**(K1+K2*np.log10(V)+K3*(np.log10(V))**2)
    cbm = cp*(B1+B2*Fm*Fp) #modulo de Bare
    c1 = cbm*(MS2015/MS2001)*f/amort #coste del reactor euros/ano
    c2 = n_OL*PM_OL*PriceOL*Horas_Ano #coste del ciclohexanol euros/ano
    c3 = n_ONA*PM_ONA*PriceONA*Horas_Ano #beneficios de la ciclohexanona euros/ano
    c4 = Vreac*PriceAceite/amort #coste del aceite termica euros/ano
    c5 = Vtub*datos.ro_l*PriceCat/amort #coste del catalizador euros/ano

    h = 1e-10

    return [c3-c1-c2-c4-c5, c1/(c1+c2+c3+c4+c5+h),c2/(c1+c2+c3+c4+c5+h),
            c3/(c1+c2+c3+c4+c5+h),c4/(c1+c2+c3+c4+c5+h),c5/(c1+c2+c3+c4+c5+h),
            c1, c2, c3, c4, c5]

Coste = np.zeros((datos.nt,datos.nl,11))

l = np.loadtxt('./data/desactivacionX.dat')
t = np.loadtxt('./data/desactivaciont.dat')+1e-6

N = np.zeros((datos.nt,datos.nl,5))
T = np.zeros((datos.nt,datos.nl,1))
P = np.zeros((datos.nt,datos.nl,1))

N[:,:,0] = np.loadtxt('./data/desactivacion0.dat') #ciclohexanol
N[:,:,1] = np.loadtxt('./data/desactivacion1.dat') #ciclohexanona
N[:,:,2] = np.loadtxt('./data/desactivacion2.dat')
N[:,:,3] = np.loadtxt('./data/desactivacion3.dat')
N[:,:,4] = np.loadtxt('./data/desactivacion4.dat')
T[:,:,0] = np.loadtxt('./data/desactivacion5.dat')
P[:,:,0] = np.loadtxt('./data/desactivacion6.dat')

for j in range(datos.nl):
    for i in range(datos.nt):
        nOL = N[0,0,0]-N[i,j,0] #el ciclo hexanol que no reacciona se reutiliza
        nONA = N[i,j,1]
        for n in range(11):
            Coste[i,j,n] = costeReactor(nOL, nONA, l[j], t[i])[n]


#representar
n0 = 1 #punto inicial

#coste total
X, Y = np.meshgrid(l, t)
fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 1, 1, projection='3d')
surf = ax.plot_surface(X[n0:,n0:], Y[n0:,n0:], Coste[n0:,n0:,0], cmap=cm.coolwarm)
plt.xlabel("longitud (m)")
plt.ylabel("tiempo (horas)")
plt.title("Coste reactor euros/ano")
plt.show()

#Relacion de costes
fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 1, 1, projection='3d')
surf = ax.plot_surface(X[n0:,n0:], Y[n0:,n0:], Coste[n0:,n0:,1], cmap=cm.coolwarm)
plt.xlabel("longitud (m)")
plt.ylabel("tiempo (horas)")
plt.title("Relacion reactor")
plt.show()

fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 1, 1, projection='3d')
surf = ax.plot_surface(X[n0:,n0:], Y[n0:,n0:], Coste[n0:,n0:,2], cmap=cm.coolwarm)
plt.xlabel("longitud (m)")
plt.ylabel("tiempo (horas)")
plt.title("Relacion ciclohaxanol")
plt.show()

fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 1, 1, projection='3d')
surf = ax.plot_surface(X[n0:,n0:], Y[n0:,n0:], Coste[n0:,n0:,3], cmap=cm.coolwarm)
plt.xlabel("longitud (m)")
plt.ylabel("tiempo (horas)")
plt.title("Relacion ciclohaxanona")
plt.show()

fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 1, 1, projection='3d')
surf = ax.plot_surface(X[n0:,n0:], Y[n0:,n0:], Coste[n0:,n0:,4], cmap=cm.coolwarm)
plt.xlabel("longitud (m)")
plt.ylabel("tiempo (horas)")
plt.title("Relacion aceite termico")
plt.show()

fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 1, 1, projection='3d')
surf = ax.plot_surface(X[n0:,n0:], Y[n0:,n0:], Coste[n0:,n0:,5], cmap=cm.coolwarm)
plt.xlabel("longitud (m)")
plt.ylabel("tiempo (horas)")
plt.title("Relacion catalizador")
plt.show()


#maximo
posicion_maximo = np.argmax(Coste[:,:,0])
posicion_maximo_2d = np.unravel_index(posicion_maximo, Coste[:,:,0].shape)
[ntmax, nLmax] = posicion_maximo_2d
print(ntmax, nLmax)
print("Longitud maxima (metros):        ",l[nLmax])
print("Tiempo maximo (horas):           ",t[ntmax])
print("Beneficio maximo (euros/ano):    ",Coste[ntmax, nLmax, 0])
print("Coste reactor (euros/ano):       ",Coste[ntmax, nLmax, 6])
print("Coste ciclohexanol (euros/ano):  ",Coste[ntmax, nLmax, 7])
print("Coste ciclohexanona (euros/ano): ",Coste[ntmax, nLmax, 8])
print("Coste aceite (euros/ano):        ",Coste[ntmax, nLmax, 9])
print("Coste catalizador (euros/ano):   ",Coste[ntmax, nLmax, 10])
