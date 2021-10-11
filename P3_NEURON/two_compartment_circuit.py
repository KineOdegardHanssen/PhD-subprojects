import math
import numpy as np 
import matplotlib.pyplot as plt 

testit = True

# Think about units... #SI
# Length: m
# Resistance: Ohm
# Voltage: V
# Loop over values of dd and ld? Or make another script that does that?

# Setting realistic paramters 
Cms = 1e-2  #Capacitance # 1mu F cm2
Cmd = 1e-2  #Capacitance # 1mu F cm2
Ra = 100e-2 # Ohm cm
Rms = 30000e-4
Rmd = 30000e-4
V0 = -65e-3 
Ie = -0.1e-9 # nA
El = -65e-3
# Soma size # In mu m
ds = 10e-6
ls = 10e-6
# Dendrite:
dd = 10e-6
ld = 10e-6

# Setting simulation parameters 
T = 1000 # Total time of simulation 
#t0 = T*100 # Length of input current
dt = 2**-11
nsteps = int(T/dt) 
L = 2 # Number of capacitors 
Vs = np.zeros(nsteps) 
Vd = np.zeros(nsteps) 
t = np.zeros(nsteps) 

daxialfacs = ds/(4*Ra*ls**2) # Need to have one for soma and one for dend
daxialfacd = dd/(4*Ra*ld**2) # Need to have one for soma and one for dend
dtdivCmds  = dt/(Cms)
Rmsinv     = 1./Rms
dtdivCmdd  = dt/(Cmd)
Rmdinv     = 1./Rmd
Iefac      = Ie/(np.pi*ds*ls)

# Starting simulation 
Vs[0] = V0 
Vd[0] = V0 
for i in range(0,nsteps-1): 
    t[i+1] = t[i] + dt 
    Vdiff  = Vd[i]-Vs[i]  # Vd[i]+Vs[i] gir "gode" resultater??
    daxials = daxialfacs*Vdiff
    daxiald = -daxialfacd*Vdiff
    dVs =((El-Vs[i])*Rmsinv+daxials+Iefac)*dtdivCmds
    Vs[i+1] = Vs[i] +dVs
    Vd[i+1] = Vd[i] +((El-Vd[i])*Rmdinv+daxiald)*dtdivCmdd
    if i<5:
        print('El-Vs[i]:',El-Vs[i], '; Vs[i]:',Vs[i])
        print('dVs:',dVs)
        print('(Vd[i]-Vs[i])*facs:',daxials)
    else:
        #Iefac=0
        if i<10:
            print('El-Vs[i]:',El-Vs[i], '; Vs[i]:',Vs[i])
            print('dVs:',dVs)
            print('(Vd[i]-Vs[i])*facs:',daxials)

print('Ie:',Ie)
print('Iefac:',Iefac)

# Plotting results 
plt.figure() 
plt.plot(t[:5],Vs[:5],label='Vsoma') 
plt.plot(t[:5],Vd[:5],label='Vdend') 
plt.xlabel('Time')
plt.ylabel('Voltage')
plt.title('Two-compartment model')
plt.legend(loc='upper right')
plt.show()

print('Vs:',Vs)
Vmax  = V0
Vmin  = min(Vs)
dV    = Vmax-Vmin
dV1e  = dV*math.exp(-1)
V1e   = Vmin+dV1e
time  = t
idelay = 0
print('Vmax:',Vmax)
print('Vmin:',Vmin)
print('V1e:',V1e)

# Will perform linear interpolation near the crossing with V1e.
tbef = 0
taft = 0
Vbef = 0
Vaft = 0
for i in range(len(Vs)):
    if Vs[i]-V1e<0:
        tbef = time[i-1]
        taft = time[i]
        Vbef = Vs[i-1]
        Vaft = Vs[i]
        print('Vs:',Vs[i-1],'i-1:',i-1,'; t:', time[i-1])
        print('Vs:',Vs[i],'i:',i,'; t:', time[i])
        break
    
a = (Vaft-Vbef)/(taft-tbef)
t_interpolated = (V1e-Vbef+a*tbef)/a-idelay
    
print('tbef',tbef)
print('taft',taft)
print('t_interpolated',t_interpolated)
    
tau1 = tbef-idelay
tau2 = taft-idelay
tau3 = t_interpolated
    
print('dV+Vmax:',dV+Vmax)
if testit==True:
    idur = T
    N = 100*idur #100000*idur #100000*idur: Very slow
    timecut = np.linspace(0,idur,N)
    fit1 = np.zeros(N)
    fit2 = np.zeros(N)
    fit3 = np.zeros(N)
    for i in range(N):
        fit1[i] = dV*math.exp(-timecut[i]/tau1)+Vmin
        fit2[i] = dV*math.exp(-timecut[i]/tau2)+Vmin
        fit3[i] = dV*math.exp(-timecut[i]/tau3)+Vmin
    timecut = np.linspace(idelay,idur+idelay,N)
    
    plt.figure(figsize=(6,5))
    plt.plot(time,Vs, label='Data')
    plt.plot(timecut,fit1, 'r--', label='tbef')
    plt.plot(timecut,fit2, 'g--', label='taft')
    plt.plot(timecut,fit3, 'b--', label='tinterpolated')
    plt.plot(tau1+idelay,V1e, 'ro', label='tbef')
    plt.plot(tau2+idelay,V1e, 'go', label='taft')
    plt.plot(tau3+idelay,V1e, 'bo', label='tinterpolated')
    plt.xlabel(r'$t$ [ms]')
    plt.ylabel(r'$V$ [mV]')
    plt.title(r'$V$ vs $t$')
    plt.tight_layout()
    plt.axis([idelay-0.5,idelay+10,min(Vs)-0.01,max(Vs)+0.01])
    plt.legend(loc='upper right')
    plt.show()

