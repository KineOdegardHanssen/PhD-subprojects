import os
from os.path import join
import sys
import matplotlib.pyplot as plt
import json
import neuron
import time as tm
import numpy as np

def give_minf_caq(v):
    N = len(v)
    minf = np.zeros(N)
    for i in range(N):
        minf[i] = 1.0/(1.0+np.exp(-(v[i]+16.3)/7.9))
    return minf

def give_minf_canin(v):
    zetam  = -3.4
    vhalfm = -21
    celsius = 37 # 34?
    alpm = np.exp(1.e-3*zetam*(v-vhalfm)*9.648e4/(8.315*(273.16+celsius)))
    minf = 1/(1+alpm)
    return minf

def give_minf_CaHVA(v):
    N = len(v)
    y = 3.8
    mInf   = np.zeros(N)
    mTau   = np.zeros(N)
    for i in range(N):
        xi = -27-v[i]
        if abs(xi/y)<1e-6:
            mAlpha = 0.055*y*(1-xi)/(2.*y) # ??? # Double check math notation in mod-files
        else:
            mAlpha = 0.055*xi/(np.exp(xi/y)-1)
        mBeta = 0.94*np.exp((-75-v[i])/17.)
        mInf[i] = mAlpha/(mAlpha+mBeta)
        mTau[i] = 1./(mAlpha+mBeta)
    return mInf, mTau

def return_Vhalf(v,minf):
    Vhalf = 0
    passedit = False
    #print('minf:',minf)
    for i in range(len(minf)):
        #print('minf[i]:',minf[i])
        if minf[i]>0.5 and passedit==False:
            Vhalf = Vs[i]
            passedit=True # Should not need this, but nice to be safe
            break
    return Vhalf

mtau = 1.13*2

Vs = np.linspace(-100,80,202)

print('Vs:',Vs)

minf_caq = give_minf_caq(Vs)
minf_canin = give_minf_canin(Vs)
minf_CaHVA, mtau_CaHVA = give_minf_CaHVA(Vs)

#print('minf_canin:',minf_canin)
print('minf_CaHVA:',minf_CaHVA)

Vhalf_caq = return_Vhalf(Vs,minf_caq)
Vhalf_canin = return_Vhalf(Vs,minf_canin)
Vhalf_CaHVA = return_Vhalf(Vs,minf_CaHVA)

plt.figure()
plt.plot(Vs,minf_caq,label=r'CaQ',color='tab:blue')
plt.plot(Vs,minf_canin,label=r'CaN',color='tab:orange')
plt.plot(Vs,minf_CaHVA,label=r'CaHVA',color='tab:green')
plt.axvline(x=Vhalf_caq,linestyle='--',linewidth=0.75,color='tab:blue')
plt.axvline(x=Vhalf_canin,linestyle='--',linewidth=0.75,color='tab:orange')
plt.axvline(x=Vhalf_CaHVA,linestyle='--',linewidth=0.75,color='tab:green')
#plt.plot(Vs,hinf,label=r'$h_\infty$')
plt.xlabel(r'$V$ (mV)')
plt.ylabel('$m_\infty$')
plt.legend(loc='lower right')
plt.title(r'$m_\infty$')
plt.show()



