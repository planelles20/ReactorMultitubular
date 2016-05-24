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
PriceONA = 1570/1000 # euros/kg
PriceOL = 1100/1000 # euros/kg
PriceAceite = 2364.52 # euros/m3
PriceCat = 200 #euros/kg
#coste de energia para calentar el refrigerante(80% de eficencia) 2001
CQ = 7.5 #$/GJ

#tiempo de amortizacion
n = 10 #tiempo aproximado que va a estar la planta en funcionamiento
i = 0.08 #interes del 8%
amort = ((1+i)**n-1)/i #anos

#horas trabajadas a la semana
horas_semana = 168
#mese al ano trabajados
meses_ano = 11
#peso molecular
PM_OL = datos.PM[0]   # kmol/kg
PM_ONA = datos.PM[1]  # kmol/kg

#Presion de diseno
P = 1*(101325/1e5)*1.10 #bares de presion

## Valores de las correlaciones
#reactor
K11 = 3.5565
K21 = 0.3776
K31 = 0.0905

Fp1 = ((P+1)*D/(2*(850-0.6*(P+1))+0.00315))/0.0063 #factor de la presion
if Fp1 < 1:
    Fp1 = 1.0

B11 = 1.49
B21 = 1.52
Fm1 = 3.1 #factor material

#intercambiador
K12 = 4.3247
K22 = -0.3030
K32 = 0.1634

C1 = 0
C2 = 0
C3 = 0

Fp2 = 10**(C1+C2*np.log10(P)+C3*(np.log10(P)**2)) #factor de la presion
if Fp2 < 1:
    Fp2 = 1.0

B12 = 1.63
B22 = 1.66
Fm2 = 2.8 #factor material

#ecuacion calcula el coste en funcion de la composicion, longitud, tiempo
# unidades: euros/ano
def costeReactor(n_OL, n_ONA, q, L, t):
    Vtub = Stub*Ntub*L #volumen de los tubos
    Vreac = Sreactor*L-Vtub #volumen que acua el aceite
    A = Cir_tub*Ntub*L #area intercambio de calor
    V = Sreactor*L+1e-10
    #horas que rabaja el rector al ano, teniendo en cuenta que necesita 8 horas
    #para regenerar  el catalizador
    if t > 8:
        Horas_Ano = meses_ano*4*(horas_semana)*((t-8)/t)
    else:
        Horas_Ano = 0

    cp1 = 10**(K11+K21*np.log10(V)+K31*(np.log10(V))**2)
    cbm1 = cp1*(B11+B21*Fm1*Fp1) #modulo de Bare para el vessel
    cp2 = 10**(K12+K22*np.log10(A)+K32*(np.log10(A))**2)
    cbm2 = cp2*(B12+B22*Fm2*Fp2) #modulo de Bare para el intercambiador
    c1 = (cbm1+cbm2)*(MS2015/MS2001)*f/amort #coste del reactor euros/ano
    c2 = n_OL*PM_OL*PriceOL*Horas_Ano #coste del ciclohexanol euros/ano
    c3 = n_ONA*PM_ONA*PriceONA*Horas_Ano #beneficios de la ciclohexanona euros/ano
    c4 = Vreac*PriceAceite/amort #coste del aceite termica euros/ano
    c5 = Vtub*datos.ro_l*PriceCat/amort #coste del catalizador euros/ano
    c6 = q*(Horas_Ano/1000/1000)*CQ**(MS2015/MS2001)*f # coste calentar refrigerante euros/ano
    h = 1e-10

    if L == 2.96593 and t == 56.000001:
        print(cp1, Fm1, Fp1, cbm1, cbm1*(MS2015/MS2001)*f)
        print(cp2, Fm2, Fp2, cbm2, cbm2*(MS2015/MS2001)*f)

    return [c3-c1-c2-c4-c5-c6, c1/(c1+c2+c3+c4+c5+c6+h),c2/(c1+c2+c3+c4+c5+c6+h),
            c3/(c1+c2+c3+c4+c5+c6+h),c4/(c1+c2+c3+c4+c5+c6+h),c5/(c1+c2+c3+c4+c5+c6+h),
            c6/(c1+c2+c3+c4+c5+c6+h), c1, c2, c3, c4, c5, c6]

Coste = np.zeros((datos.nt,datos.nl,13))

l = np.loadtxt('./data/desactivacionX.dat')
t = np.loadtxt('./data/desactivaciont.dat')+1e-6

N  = np.zeros((datos.nt,datos.nl,5))
T  = np.zeros((datos.nt,datos.nl,1))
Ts = np.zeros((datos.nt,datos.nl,1))
P  = np.zeros((datos.nt,datos.nl,1))
Q  = np.zeros((datos.nt,datos.nl,1)) #calor kJ/h

N[:,:,0]  = np.loadtxt('./data/desactivacion0.dat') #ciclohexanol
N[:,:,1]  = np.loadtxt('./data/desactivacion1.dat') #ciclohexanona
N[:,:,2]  = np.loadtxt('./data/desactivacion2.dat')
N[:,:,3]  = np.loadtxt('./data/desactivacion3.dat')
N[:,:,4]  = np.loadtxt('./data/desactivacion4.dat')
T[:,:,0]  = np.loadtxt('./data/desactivacion5.dat')
P[:,:,0]  = np.loadtxt('./data/desactivacion7.dat')
Ts[:,:,0] = np.loadtxt('./data/desactivacion6.dat')

#calculo del calor en cada punto
for j in range(datos.nl):
    for i in range(datos.nt):
        if j == 0:
            Q[i,j] = datos.U*(Ts[i,j]-T[i,j])*(Cir_tub*datos.Ntub*(l[1]-l[0]))
        else:
            Q[i,j] = datos.U*(Ts[i,j]-T[i,j])*(Cir_tub*datos.Ntub*(l[1]-l[0])) + Q[i,j-1]

#calculo del los beneficios, costes...
for j in range(datos.nl):
    for i in range(datos.nt):
        nOL = N[0,0,0]#-N[i,j,0] #el ciclo hexanol que no reacciona se reutiliza
        nONA = N[i,j,1]
        for nn in range(13):
            q = Q[i,j] #calor en el pinto i,j en kJ/h
            Coste[i,j,nn] = costeReactor(nOL, nONA, q, l[j], t[i])[nn]



#maximo
posicion_maximo = np.argmax(Coste[:,:,0])
posicion_maximo_2d = np.unravel_index(posicion_maximo, Coste[:,:,0].shape)
[ntmax, nLmax] = posicion_maximo_2d
print(ntmax, nLmax)
print("Longitud maxima (metros):        ",l[nLmax])
print("Tiempo maximo (horas):           ",t[ntmax])
print("Beneficio maximo (euros/ano):    ",Coste[ntmax, nLmax, 0])
print("Coste reactor (euros):           ",Coste[ntmax, nLmax, 7]*amort)
print("Coste ciclohexanol (euros/ano):  ",Coste[ntmax, nLmax, 8])
print("Ventas ciclohexanona (euros/ano):",Coste[ntmax, nLmax, 9])
print("Coste aceite euros:              ",Coste[ntmax, nLmax, 10]*amort)
print("Coste catalizador euros:         ",Coste[ntmax, nLmax, 11]*amort)
print("Coste calor (euros/ano):         ",Coste[ntmax, nLmax, 12])
print("Calor (kJ/h):                    ",Q[ntmax,nLmax,0])

#representar
n0 = 1 #punto inicial

#coste total
X, Y = np.meshgrid(l, t)
fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 1, 1, projection='3d')
surf = ax.plot_surface(X[n0:,n0:], Y[n0:,n0:], Coste[n0:,n0:,0], cmap=cm.coolwarm)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.xlabel("longitud (m)")
plt.ylabel("tiempo (horas)")
plt.title("Beneficios reactor euros/ano")
plt.show()

fig = plt.figure(1)
maximo = Coste[ntmax, nLmax, 0]
levels = [-6e7,-5e7,-4e7,-3e7,-2e7,-1e7, -1e6, -1e5, -1e4, -1e3, -1e2, -1e1, 0, maximo/10, maximo/9, maximo/8, maximo/7, maximo/6, maximo/5, maximo/4, maximo/3, maximo/2, maximo/1.5, maximo/1.15, maximo/1.10, maximo/1.08, maximo/1.05, maximo/1.03, maximo/1.01, maximo/1.001, maximo/1.0001, maximo/1.00001, maximo]

plt.contour(X[n0:,n0:], Y[n0:,n0:], Coste[n0:,n0:,0], levels, cmap=cm.coolwarm)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.xlabel("longitud (m)")
plt.ylabel("tiempo (horas)")
plt.title("Beneficios reactor euros/ano")
plt.show()


fig = plt.figure(1)
plt.title("")
plt.xlabel("longitud (m)")
plt.ylabel("Beneficios reactor €/año")
plt.plot(l[n0:], Coste[ntmax,n0:,0])
plt.show()

fig = plt.figure(1)
plt.title("")
plt.xlabel("tiempo (h)")
plt.ylabel("Beneficios reactor €/año")
plt.plot(t[n0:], Coste[n0:,nLmax,0])
plt.show()

#estudio economico

##              P-100     E-100/E-101  E-102    E-103      V-100       reactor
CosteEquipos = [36356.60, 716199.52, 266946.82, 125220.84, 65458.73, Coste[ntmax, nLmax, 7]*amort]
#                 P-100   E-102   E-103   reactor
CosteOperacion = [1.00e4, 5.42e4, 1.00e4, Coste[ntmax, nLmax, 12]]

#coste total planta
CTM = 4.74*sum(CosteEquipos)
# coste
CGR = CTM+0.5*sum(CosteEquipos)

def CostFlow(t):
    #coste capital total
    FCI = Coste[ntmax, nLmax, 11]*amort+Coste[ntmax, nLmax, 10]*amort+CGR
    #Net earnings
    NE = 0
    if t > 0:
        Ingresos = Coste[ntmax, nLmax, 9]
        Gastos =  Coste[ntmax, nLmax, 8]+sum(CosteOperacion)
        Impuestos = 0.2*(Ingresos-Gastos-0.15*FCI)
        NE = Ingresos-Gastos-Impuestos
    #devolucion
    DEV = 0
    if t == n:
        DEV = 0.1*FCI

    if t > 0:
        FCI = 0

    CF = NE-FCI+DEV

    return CF

y = np.zeros((n+1))

for t in range(n+1):
    if t == 0:
        y[t] = CostFlow(t)
    else:
        y[t] = CostFlow(t)+y[t-1]

P = np.zeros((n+2))
for t in range(n+2):
    CF = 0
    if t >= 2:
        #coste capital total
        FCI = Coste[ntmax, nLmax, 11]*amort+Coste[ntmax, nLmax, 10]*amort+CGR
        #Net earnings
        NE = 0
        if t > 0:
            Ingresos = Coste[ntmax, nLmax, 9]
            Gastos =  Coste[ntmax, nLmax, 8]+sum(CosteOperacion)
            Impuestos = 0.2*(Ingresos-Gastos-0.15*FCI)
            NE = Ingresos-Gastos-Impuestos
        #devolucion
        DEV = 0
        if t == n+1:
            DEV = 0.1*FCI

        if t > 0:
            FCI = 0

        CF = NE+DEV

    P[t] = CF/((1+0.04)**t)

print("VAN (€): ", sum(P)+y[0])

fig = plt.figure(1)
x = [0, 2, 3, 4, 5, 6, 7, 8 ,9, 10, 11, 12]
yy = np.zeros((n+2))
for i in range(n+2):
    if i >= 1:
        yy[i] = y[i-1]

plt.title("Coste acumulado")
plt.xlabel("tiempo (años)")
plt.ylabel("Coste anualizado (€/año)")
plt.plot(x, yy)
plt.plot(x, np.ones((n+2)))
plt.show()
