import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, quad
import os
from math import *
import scipy.interpolate
from string import *
from  pylab  import *
from  mpl_toolkits.mplot3d  import Axes3D


## data

x0 = 10000
c = 1.2
t = np.linspace(0,100,800)
K = 30000
q = 0.2
p = 15000
k = 50000

## 1 et 2 : fishing = h

lh = [0,1000, 10000, 11990, 12000, 12010, 13000]

for h in lh :
    y = (x0 - h/c) * np.exp(c*t) + h/c
    plt.plot(t,y,label= h)
plt.xlim(0,15) # zoom
plt.ylim(0,2000000)
plt.title('Biomass evolution, depending of time and parameter h')
plt.xlabel('time (in years)')
plt.ylabel('biomass (en tons)')
plt.grid()
plt.legend()
plt.show()


## 3: Plot of the integral curve of Verhulst

x0 = 10000
y = x0/((1-x0/K)*(np.exp(-c*t)+(x0/(K-x0))))
plt.plot(t,y, label = "biomass with  x0 = 10000")
x0 = 50000
z = x0/((1-x0/K)*(np.exp(-c*t)+(x0/(K-x0))))
plt.plot(t,z, label = "biomass with c x0 = 50000")
plt.xlim(0,7) # zoom
plt.ylim(9000,51000)
plt.title('Biomass evolution, depending of time')
plt.xlabel('time (in years)')
plt.ylabel('biomass (in tons)')
plt.grid()
plt.legend()
plt.show()

## 4 Plot of the integral curve of Verhulst with h different of 0
lh = [0,c*K/8,7990, c*K/4 - 10, c*K/4, c*K/4 + 10, c*K/3]
x0 = 10000
for h in lh :
    def g(x,t):
        return c*x*(1-x/K) - h
    j =  odeint(g, x0, t)
    plt.plot(t, j, label = (h, x0))

y0 = 50000
for h in lh :
    def g(x,t):
        return c*x*(1-x/K) - h
    j =  odeint(g, y0, t)
    plt.plot(t, j, label = (h, y0))
plt.xlim(0,30) # zoom
plt.ylim(0,51000)
plt.title('Biomass evolution, depending of time')
plt.xlabel('time (in years)')
plt.ylabel('biomass (in tons)')
plt.grid()
plt.legend()
plt.show()

## 6 Plot of the integral curve of Verhulst with h different of 0 and y0 = 16 000
lh = [0,c*K/8,7990, c*K/4 - 10, c*K/4, c*K/4 + 10, c*K/3]
x0 = 10000
for h in lh :
    def g(x,t):
        return c*x*(1-x/K) - h
    j =  odeint(g, x0, t)
    plt.plot(t, j, label = (h, x0))

y0 = 16000
for h in lh :
    def g(x,t):
        return c*x*(1-x/K) - h
    j =  odeint(g, y0, t)
    plt.plot(t, j, label = (h, y0))
plt.xlim(0,10) # zoom
plt.ylim(0,31000)
plt.title('Biomass evolution, depending of time')
plt.xlabel('time (in years)')
plt.ylabel('biomass (in tons)')
plt.grid()
plt.legend()
plt.show()

## 5: Plot of the MSY curve

x = np.linspace(0, 30000, 30000)
plt.plot(x,x*c*(1-x/K),label= "Stock growth")
plt.xlim(0,30001) # zoom
plt.ylim(0,10000)
plt.title('Growing depending of the stock size')
plt.xlabel('Biomass (en tons)')
plt.ylabel('Growing of the biomass (in ton per year)')
plt.grid()
plt.legend()
plt.show()


## 7 Plot of the 3D curve with qEx

x0 = 20000
ax = Axes3D( figure () )
X = np.arange (0,  50 ,  0.5)
Y = np.arange (0,  7,  0.5)
t,E = np.meshgrid (X, Y)
Z = (x0*(1-q*E/c)*np.exp(t*(c-q*E)))/(1-q*E/c + x0/K *(np.exp(t*(c-q*E))-1))
ax.plot_surface (t, E,  Z , cmap='cool', label = "x0 = 20 000" )
x0 = 10000
Z = (x0*(1-q*E/c)*np.exp(t*(c-q*E)))/(1-q*E/c + x0/K *(np.exp(t*(c-q*E))-1))
ax.plot_surface (t, E,  Z , cmap='hot', label = "x0 = 10 000" )
plt.title("Biomass evolution, depending of time and fishing effort for an exploitated specie (Verhulst with h(t,x) = qEx) (Verhulst avec h(t,x) = qEx)")
ax.autoscale_view(tight=None, scalex=True, scaley=True, scalez=True)
plt.xlabel('time (in years)')
plt.ylabel('Fishing effort (per year)')
show()


## 8 : Income depending of time and effort
x0 = 20000
ax = Axes3D( figure () )
X = np.arange (0,  50 ,  0.5)
Y = np.arange (0,  7,  0.5)
t,E = np.meshgrid(X, Y)
Z = E*(p*q*((x0*(1-q*E/c)*np.exp(t*(c-q*E)))/(1-q*E/c + x0/K *(np.exp(t*(c-q*E))-1))) - k)
plt.xlabel('time (in years)')
plt.ylabel('Fishing effort (per year)')
plt.title("Profits (in euros) depending of time and fishing effort")
ax.autoscale_view(tight= None, scalex=True, scaley=True, scalez=True)
ax.plot_surface (t, E,  Z , cmap='cool' )
show()

## 9 Ploting Geq(x)

x = np.linspace(0, 30000, 1000)
plt.plot(x, c*(p*x- k/q)*(1-x/k))
plt.title("Profits (in euros) at equilibrium depending of x")
plt.xlabel('biomass (en tons)')
plt.ylabel('Profit (in euros per year)')
plt.grid()
plt.legend()
plt.show()

## 10 Ploting x delta

x = np.linspace(0, 2, 1000)
plt.plot(x, (K/(4*p*q))* (k/K + p*q*(1-x) + np.sqrt((k/K + p*q*(1-x))**2 + 8*p*q*k*x/K)))
plt.xlabel('δ/c')
plt.ylabel('x(δ/c) (in tons)')
plt.grid()
plt.legend()
plt.show()
