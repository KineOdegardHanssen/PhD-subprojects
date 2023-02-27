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
m_standard, mInf_standard, mTau = mInf_graph(vs,xstandard,ystandard)
m_st_p5, mInf_st_p5, mTau_p5 = mInf_graph(vs+5,xstandard-5,ystandard)
m_st_p10, mInf_st_p10, mTau_p10 = mInf_graph(vs+10,xstandard-10,ystandard)
m_st_p50, mInf_st_p50, mTau_p50 = mInf_graph(vs+50,xstandard-50,ystandard)
h_standard, hInf_standard, hTau = hInf_graph(vs)
relcond_standard = relcond(m_standard,h_standard)

print('max(mTau):',max(mTau),'; min(mTau):',min(mTau))
print('max(hTau):',max(hTau),'; min(hTau):',min(hTau))


xvp10 = -17-vs
xvm10 = -37-vs
yp1 = 4.8
ym1 = 2.8
m_vp10, mInf_vp10, mTau = mInf_graph(vs,xvp10,ystandard)
m_vm10, mInf_vm10, mTau = mInf_graph(vs,xvm10,ystandard)
m_yp1, mInf_yp1, mTau = mInf_graph(vs,xstandard,yp1)
m_ym1, mInf_ym1, mTau = mInf_graph(vs,xstandard,ym1)
relcond_vp10 = relcond(m_vp10,h_standard)
relcond_vm10 = relcond(m_vm10,h_standard)
relcond_yp1  = relcond(m_yp1,h_standard)
relcond_ym1  = relcond(m_ym1,h_standard)

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

m_st_plain_p5, mInf_st_plain_p5, mTau_plain_p5 = mInf_graph_plain(vs+5,ystandard)
m_st_plain_p10, mInf_st_plain_p10, mTau_plain_p10 = mInf_graph_plain(vs+10,ystandard)
m_st_plain_p50, mInf_st_plain_p50, mTau_plain_p50 = mInf_graph_plain(vs+50,ystandard)
m_st_plain_m5, mInf_st_plain_m5, mTau_plain_m5 = mInf_graph_plain(vs-5,ystandard)
m_st_plain_m10, mInf_st_plain_m10, mTau_plain_m10 = mInf_graph_plain(vs-10,ystandard)
m_st_plain_m50, mInf_st_plain_m50, mTau_plain_m50 = mInf_graph_plain(vs-50,ystandard)

plt.figure(figsize=(6,5))
plt.plot(vs,mInf_standard,color='k',label='mInf(V)')
plt.plot(vs,mInf_st_p5,color='r',label='mInf(V+5mV)')
plt.plot(vs,mInf_st_p10,color='b',label='mInf(V+10mV)')
plt.plot(vs,mInf_st_p50,color='g',label='mInf(V+50mV)')
#plt.plot(vs,hInf_standard,color='tab:red',label='hInf')
###plt.plot(vs,relcond_vm10,color='tab:pink',label='vm10 mAlpha')
###plt.plot(vs,relcond_yp1,color='darkblue',label='yp1 mAlpha')
###plt.plot(vs,relcond_ym1,color='tab:blue',label='ym1 mAlpha')
plt.xlabel('V [mV]')
plt.ylabel(r'Activation')
plt.legend(loc='center right')
plt.tight_layout()
plt.show()

plt.figure(figsize=(6,5))
plt.plot(vs,mInf_standard,color='k',label='mInf(V)')
plt.plot(vs,mInf_st_plain_p5,color='r',label='mInf(V+5mV)')
plt.plot(vs,mInf_st_plain_p10,color='b',label='mInf(V+10mV)')
plt.plot(vs,mInf_st_plain_p50,color='g',label='mInf(V+50mV)')
plt.plot(vs,mInf_st_plain_m5,color='m',label='mInf(V-5mV)')
plt.plot(vs,mInf_st_plain_m10,color='c',label='mInf(V-10mV)')
plt.plot(vs,mInf_st_plain_m50,color='lightgreen',label='mInf(V-50mV)')
#plt.plot(vs,hInf_standard,color='tab:red',label='hInf')
###plt.plot(vs,relcond_vm10,color='tab:pink',label='vm10 mAlpha')
###plt.plot(vs,relcond_yp1,color='darkblue',label='yp1 mAlpha')
###plt.plot(vs,relcond_ym1,color='tab:blue',label='ym1 mAlpha')
plt.xlabel('V [mV]')
plt.ylabel(r'Activation')
plt.legend(loc='center right')
plt.title('Plain function')
plt.tight_layout()
plt.show()


vs = np.linspace(-80,40,1000)

m_standard, mInf_standard, mTau = mInf_graph(vs,xstandard,ystandard)
m_st_p5, mInf_st_p5, mTau_p5 = mInf_graph(vs+5,xstandard-5,ystandard)
m_st_p10, mInf_st_p10, mTau_p10 = mInf_graph(vs+10,xstandard-10,ystandard)
m_st_p50, mInf_st_p50, mTau_p50 = mInf_graph(vs+50,xstandard-50,ystandard)
h_standard, hInf_standard, hTau = hInf_graph(vs)

m_st_plain_p5, mInf_st_plain_p5, mTau_plain_p5 = mInf_graph_plain(vs+5,ystandard)
m_st_plain_p10, mInf_st_plain_p10, mTau_plain_p10 = mInf_graph_plain(vs+10,ystandard)
m_st_plain_p50, mInf_st_plain_p50, mTau_plain_p50 = mInf_graph_plain(vs+50,ystandard)
m_st_plain_m5, mInf_st_plain_m5, mTau_plain_m5 = mInf_graph_plain(vs-5,ystandard)
m_st_plain_m10, mInf_st_plain_m10, mTau_plain_m10 = mInf_graph_plain(vs-10,ystandard)
m_st_plain_m50, mInf_st_plain_m50, mTau_plain_m50 = mInf_graph_plain(vs-50,ystandard)

plt.figure(figsize=(6,5))
plt.plot(vs,mInf_standard,color='k',label='mInf')
plt.plot(vs,mInf_st_p5,color='r',label='mInf(V+5mV)')
plt.plot(vs,mInf_st_p10,color='b',label='mInf(V+10mV)')
plt.plot(vs,mInf_st_p50,color='g',label='mInf(V+50mV)')
#plt.plot(vs,hInf_standard,color='tab:red',label='hInf')
##plt.plot(vs,relcond_vm10,color='tab:pink',label='vm10 mAlpha')
###plt.plot(vs,relcond_yp1,color='darkblue',label='yp1 mAlpha')
###plt.plot(vs,relcond_ym1,color='tab:blue',label='ym1 mAlpha')
plt.xlabel('V [mV]')
plt.ylabel(r'Activation')
plt.legend(loc='center right')
plt.tight_layout()
plt.show()


plt.figure(figsize=(6,5))
plt.plot(vs,mInf_standard,color='k',label='mInf(V)')
plt.plot(vs,mInf_st_plain_p5,color='r',label='mInf(V+5mV)')
plt.plot(vs,mInf_st_plain_p10,color='b',label='mInf(V+10mV)')
plt.plot(vs,mInf_st_plain_p50,color='g',label='mInf(V+50mV)')
plt.plot(vs,mInf_st_plain_m5,color='m',label='mInf(V-5mV)')
plt.plot(vs,mInf_st_plain_m10,color='c',label='mInf(V-10mV)')
plt.plot(vs,mInf_st_plain_m50,color='lightgreen',label='mInf(V-50mV)')
#plt.plot(vs,hInf_standard,color='tab:red',label='hInf')
###plt.plot(vs,relcond_vm10,color='tab:pink',label='vm10 mAlpha')
###plt.plot(vs,relcond_yp1,color='darkblue',label='yp1 mAlpha')
###plt.plot(vs,relcond_ym1,color='tab:blue',label='ym1 mAlpha')
plt.xlabel('V [mV]')
plt.ylabel(r'Activation')
plt.legend(loc='center right')
plt.title('Plain function')
plt.tight_layout()
plt.show()
