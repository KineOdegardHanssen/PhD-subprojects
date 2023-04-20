import os
from os.path import join
import sys
import matplotlib.pyplot as plt
import json
import neuron
import time as tm
import numpy as np

def give_minf(v):
    minf = 1/(1+np.exp((v-(-16.3))/(-7.9)))
    return minf

mtau = 1.13*2

Vs = np.linspace(-100,80,202)

minf = give_minf(Vs)

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
plt.savefig('all_minfty_together.png')
