import matplotlib.pyplot as plt
import numpy as np

def mInf_graph(v,x,y):
    N = len(v)
    m0 = 0 # Not so sure about this one
    mInf   = np.zeros(N)
    mTau   = np.zeros(N)
    m      = np.zeros(N)
    for i in range(N):
        xi = x[i]
        if abs(xi/y)<1e-6:
            mAlpha = 0.055*y*(1-xi)/(2.*y) # ??? # Double check math notation in mod-files
        else:
            mAlpha = 0.055*xi/(np.exp(xi/y)-1)
        mBeta = 0.94*np.exp((-75-v[i])/17.)
        mInf[i] = mAlpha/(mAlpha+mBeta)
        mTau[i] = 1./(mAlpha+mBeta)
        if i==0:
            m[i]    = mInf[i]-(mInf[i]-m0)*np.exp(-300/mTau[i])
        else:
            m[i]    = mInf[i]-(mInf[i]-m[i-1])*np.exp(-300/mTau[i])
    return m, mInf, mTau

def hInf_graph(v):
    N = len(v)
    h0 = 0 # Really not sure about this one either...
    hInf   = np.zeros(N)
    hTau   = np.zeros(N)
    h      = np.zeros(N)
    for i in range(N):
        vi = v[i]
        hAlpha = 0.000457*np.exp(-(13+vi)/50.)
        hBeta = 0.0065/(np.exp(-(15+vi)/28.)+1)
        hInf[i] = hAlpha/(hAlpha+hBeta)
        hTau[i] = 1./(hAlpha+hBeta)
        if i==0:
            h[i]    = hInf[i]-(hInf[i]-h0)*np.exp(-300/hTau[i])
        else:
            h[i]    = hInf[i]-(hInf[i]-h[i-1])*np.exp(-300/hTau[i])
    return h, hInf, hTau

def zInf_graph(cai):
    N = len(cai)
    zInf = np.zeros(N)
    zInf = 1/(1+(0.00043/cai)**4.8)
    return zInf

cais = np.logspace(-5, -2, num=100, endpoint=True, base=10.0, dtype=None, axis=0)
zInf_standard = zInf_graph(cais)


caihalf_SK = 0
passedit_SK = False
for i in range(len(zInf_standard)):
    if zInf_standard[i]>0.5 and passedit_SK==False:
        caihalf_SK = cais[i]
        passedit_SK=True # Should not need this, but nice to be safe
        break

vs = np.linspace(-150,150,1000)

xstandard = -27-vs
ystandard = 3.8
m_standard, mInf_standard, mTau_standard = mInf_graph(vs,xstandard,ystandard)
h_standard, hInf_standard, hTau_standard = hInf_graph(vs)

Vhalf = 0
passedit = False
for i in range(len(mInf_standard)):
    if mInf_standard[i]>0.5 and passedit==False:
        Vhalf = vs[i]
        passedit=True # Should not need this, but nice to be safe
        break

print('max(mTau_standard):',max(mTau_standard),'; min(mTau_standard):',min(mTau_standard))
print('max(hTau_standard):',max(hTau_standard),'; min(hTau_standard):',min(hTau_standard))

print('caihalf_SK:',caihalf_SK)
plt.figure(figsize=(5,2.7),dpi=300)
plt.plot(cais,zInf_standard,color='k')
plt.axvline(x=1e-4,color='grey',linestyle=':',linewidth=0.75)
plt.axvline(x=caihalf_SK,color='k',linestyle='--',linewidth=0.75)
plt.xscale("log")
plt.xlabel(r'[Ca$^{2+}$]$_\mathregular{in}$ (mM)',fontsize=12)
plt.ylabel(r'Activation',fontsize=12)
plt.title(r'$z_\infty$',fontsize=18)
plt.tight_layout()
plt.savefig('zInf_SK.png')
plt.show()

'''
plt.figure(figsize=(6,5))
plt.plot(vs,relcond_standard,color='k',label='standard mAlpha')
plt.plot(vs,relcond_vp10,color='tab:red',label='vp10 mAlpha')
plt.plot(vs,relcond_vm10,color='tab:pink',label='vm10 mAlpha')
plt.plot(vs,relcond_yp1,color='darkblue',label='yp1 mAlpha')
plt.plot(vs,relcond_ym1,color='tab:blue',label='ym1 mAlpha')
plt.xlabel('V [mV]')
plt.ylabel(r'$\bar{g}/\bar{g}_{\mathregular{max}}$')
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()
'''

'''

fig, (ax1, ax2) = plt.subplots(2,1,figsize=(5,7),dpi=300)
ax1.plot(vs,mInf_standard,color='k')
ax1.axvline(x=Vhalf,color='k',linestyle='--',linewidth=0.75)
#ax1.set_xlabel('V (mV)')
ax1.set_ylabel(r'Activation')
ax1.set_title(r'$m_\infty$',fontsize=12)
ax1.set_title('A', loc='left',fontsize=16)

ax2.plot(vs,mTau_standard,color='k')
ax2.axvline(x=Vhalf,color='k',linestyle='--',linewidth=0.75)
ax2.set_xlabel('V (mV)')
ax2.set_ylabel(r'Time constant (ms)')
ax2.set_title(r'$\tau_m$',fontsize=12)
ax2.set_title('B', loc='left',fontsize=16)

plt.tight_layout()
plt.savefig('activation_minf_mtau.png')
plt.show()
'''
