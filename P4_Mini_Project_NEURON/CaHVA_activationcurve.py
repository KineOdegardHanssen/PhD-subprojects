import matplotlib.pyplot as plt
import numpy as np

# ... Should return mInf and hInf instead...
def mAlpha_graph(x,y):
    N = len(x)
    mAlpha = np.zeros(N)
    for i in range(N):
        xi = x[i]
        if abs(xi/y)<1e-6:
            mAlpha[i] = 0.055*y*(1-xi)/(2.*y) # ??? # Double check math notation in mod-files
        else:
            mAlpha[i] = 0.055*xi/(np.exp(xi/y)-1)
    return mAlpha

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

def mInf_graph_plain(v,y):
    N = len(v)
    m0 = 0 # Not so sure about this one
    mInf   = np.zeros(N)
    mTau   = np.zeros(N)
    m      = np.zeros(N)
    for i in range(N):
        xi = -27-v[i]
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

def relcond(m,h):
    N = len(m)
    cond = np.zeros(N)
    for i in range(N):
        cond[i] = m[i]*m[i]*h[i]
    return cond

vs = np.linspace(-150,150,1000)

xstandard = -27-vs
ystandard = 3.8
malpha_standard = mAlpha_graph(xstandard,ystandard)
m_standard, mInf_standard, mTau_standard = mInf_graph(vs,xstandard,ystandard)
h_standard, hInf_standard, hTau_standard = hInf_graph(vs)
relcond_standard = relcond(m_standard,h_standard)


m_standard, mInf_standard, mTau_standard = mInf_graph_plain(vs,ystandard)

Vhalf = 0
passedit = False
for i in range(len(mInf_standard)):
    if mInf_standard[i]>0.5 and passedit==False:
        Vhalf = vs[i]
        passedit=True # Should not need this, but nice to be safe
        break

print('max(mTau_standard):',max(mTau_standard),'; min(mTau_standard):',min(mTau_standard))
print('max(hTau_standard):',max(hTau_standard),'; min(hTau_standard):',min(hTau_standard))


xvp10 = -17-vs
xvm10 = -37-vs
yp1 = 4.8
ym1 = 2.8
vs_shifted = vs-14.5
x_shifted  = xstandard-14.5
m_vp10, mInf_vp10, mTau_vp10 = mInf_graph(vs,xvp10,ystandard)
m_vm10, mInf_vm10, mTau_vm10 = mInf_graph(vs,xvm10,ystandard)
m_shifted, mInf_shifted, mTau_shifted = mInf_graph(vs_shifted,x_shifted,ystandard)
m_yp1, mInf_yp1, mTau_yp1 = mInf_graph(vs,xstandard,yp1)
m_ym1, mInf_ym1, mTau_ym1 = mInf_graph(vs,xstandard,ym1)
relcond_vp10 = relcond(m_vp10,h_standard)
relcond_vm10 = relcond(m_vm10,h_standard)
relcond_yp1  = relcond(m_yp1,h_standard)
relcond_ym1  = relcond(m_ym1,h_standard)


m_st_plain_p14p5, mInf_st_plain_p14p5, mTau_plain_p14p5 = mInf_graph_plain(vs-14.5,ystandard)


Vhalf_shifted = 0
passedit = False
for i in range(len(mInf_standard)):
    if mInf_st_plain_p14p5[i]>0.5 and passedit==False:
        Vhalf_shifted = vs[i]
        passedit=True # Should not need this, but nice to be safe
        break

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
print('Vhalf:',Vhalf)

plt.figure(figsize=(6,5))
plt.plot(vs,mInf_standard,color='k',label='mInf')
plt.plot(vs,hInf_standard,color='tab:red',label='hInf')
#plt.plot(vs,relcond_vm10,color='tab:pink',label='vm10 mAlpha')
#plt.plot(vs,relcond_yp1,color='darkblue',label='yp1 mAlpha')
#plt.plot(vs,relcond_ym1,color='tab:blue',label='ym1 mAlpha')
plt.xlabel('V [mV]')
plt.ylabel(r'Activation')
plt.legend(loc='center right')
plt.tight_layout()
plt.show()

fig, (ax1, ax2) = plt.subplots(2,1,figsize=(5,5),dpi=300)
ax1.plot(vs,mInf_standard,color='k',label=r'$V_\mathregular{shift}=0$ mV')
ax1.plot(vs,mInf_st_plain_p14p5,color='tab:blue',label=r'$V_\mathregular{shift}=14.5$ mV')
ax1.axvline(x=Vhalf,color='k',linestyle='--',linewidth=0.75)
ax1.axvline(x=Vhalf_shifted,color='tab:blue',linestyle='--',linewidth=0.75)
#ax1.set_xlabel('V (mV)')
ax1.set_ylabel(r'Activation',fontsize=12)
ax1.set_title(r'$m_\infty$',fontsize=16)
ax1.set_title('A', loc='left',fontsize=18)
ax1.legend(loc='lower right')

ax2.plot(vs,mTau_standard,color='k')
ax2.axvline(x=Vhalf,color='k',linestyle='--',linewidth=0.75)
ax2.axvline(x=Vhalf_shifted,color='tab:blue',linestyle='--',linewidth=0.75)
ax2.set_xlabel('V (mV)',fontsize=12)
ax2.set_ylabel(r'Time constant (ms)',fontsize=12)
ax2.set_title(r'$\tau_m$',fontsize=16)
ax2.set_title('B', loc='left',fontsize=18)

plt.tight_layout()
plt.savefig('CaHVA_activation_minf_mtau.png')
plt.show()

vs = np.linspace(-80,40,1000)

m_standard, mInf_standard, mTau = mInf_graph(vs,xstandard,ystandard)
h_standard, hInf_standard, hTau = hInf_graph(vs)

plt.figure(figsize=(6,5))
plt.plot(vs,mInf_standard,color='k',label='mInf')
plt.plot(vs,hInf_standard,color='tab:red',label='hInf')
#plt.plot(vs,relcond_vm10,color='tab:pink',label='vm10 mAlpha')
#plt.plot(vs,relcond_yp1,color='darkblue',label='yp1 mAlpha')
#plt.plot(vs,relcond_ym1,color='tab:blue',label='ym1 mAlpha')
plt.xlabel('V [mV]')
plt.ylabel(r'Activation')
plt.legend(loc='center right')
plt.tight_layout()
plt.show()
