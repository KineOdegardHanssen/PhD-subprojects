import os
from os.path import join
import sys
import matplotlib.pyplot as plt
import json
import neuron
import time as tm
import numpy as np



def give_minf(v):
    alpm = np.exp(1.e-3*zetam*(v-vhalfm)*9.648e4/(8.315*(273.16+celsius)))
    minf = 1/(1+alpm)
    return minf


def give_hinf(v):
    alph = np.exp(1.e-3*zetah*(v-vhalfh)*9.648e4/(8.315*(273.16+celsius))) 
    hinf = 1/(1+alph)
    return hinf


zetam  = -3.4
zetah  = 2
vhalfm = -21
vhalfh = -40
celsius = 37 # 34?

Vs = np.linspace(-100,80,202)

minf = give_minf(Vs)
hinf = give_hinf(Vs)

Vhalf = 0
passedit = False
for i in range(len(minf)):
    if minf[i]>0.5 and passedit==False:
        Vhalf = Vs[i]
        passedit=True # Should not need this, but nice to be safe
        break


print('Vhalf:',Vhalf)

plt.figure()
plt.plot(Vs,minf,color='k',label=r'$m_\infty$')
plt.axvline(x=Vhalf,color='k',linestyle='--',linewidth=0.75)
#plt.plot(Vs,hinf,label=r'$h_\infty$')
plt.xlabel(r'$V$ (mV)')
plt.ylabel('Activation')
plt.legend(loc='center right')
plt.title(r'$m_\infty$')
plt.show()





