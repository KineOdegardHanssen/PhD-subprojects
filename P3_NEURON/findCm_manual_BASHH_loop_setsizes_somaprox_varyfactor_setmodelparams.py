import numpy as np
import matplotlib.pyplot as plt
import math

dtexp = -8

testit = False # True # # If I am unsure about the results, I can check the fit.

varymech = 'Na' # 'K' # 'leak'
varyE_bool = True
varyE = 30 #[30,40,50,60,70]
varyg = 'None'
    
varylist = [] # Should be redundant
plotstring = '_vary'
if varyE_bool==True:
    #varylist = varyE
    plotstring = plotstring + 'E'
else:
    #varylist = varyg
    plotstring = plotstring + 'g'
Nvary    = len(varylist)

if varymech=='Na':
    varyfolder = 'VaryNa/'
elif varymech=='K':
    varyfolder = 'VaryK/'
elif varymech=='leak':
    varyfolder = 'VaryLeak/'

changestring =''
if varyE_bool==True:
    changestring = changestring+'_E'+str(varyE)+'_gdflt'
else:
    changestring = changestring+'_Edefault_g'+str(varyg)

# loop de loop? # I really should...
somasize = 10
dendlen  = 1000
denddiam = 2 
L      = somasize
diam   = somasize
A      = np.pi*L*diam
Ra     = 100
v_init = -65
gpas   = 0.0003
vpas   = -65
long   = True
    
# Values of membrane capacitance (parameters):
cmfacs = [0.1,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2,3,5,6,6.5,7,8,8.8,10]#[0.1,0.5,0.75,0.88,0.94,1.0,1.25,1.5,2.0,3.0]
outcms  = []
outtaus = []
Rins    = []

# Change current:
idur = 1000   # ms
iamp = -0.1  # nA 
idelay = 100

folder = 'Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+varyfolder +'current_idur%i_iamp'% idur+str(iamp)+'/' 
outfilename = folder +'basHH_cmfs_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+changestring+'_Ctot'
outfilename_tau = folder +'basHH_cmfs_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+changestring+'_tau'
outfilename_R = folder +'basHH_cmfs_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+changestring+'_Rin'
if long==True:
    plotname        = outfilename + '_highCms.png'
    plotname_tau    = outfilename_tau + '_highCms.png'
    plotname_R      = outfilename_R + '_highCms.png'
    outfilename     = outfilename + '_highCms.txt'
    outfilename_tau = outfilename_tau + '_highCms.txt'
else:
    plotname        = outfilename + '.png'
    plotname_tau    = outfilename_tau + '.png'
    plotname_R      = outfilename_R + '.png'
    outfilename     = outfilename + '.txt'
    outfilename_tau = outfilename_tau + '.txt'
outfile = open(outfilename,'w')
outfile_R = open(outfilename_R,'w')
outfile_tau = open(outfilename_tau,'w')

for cmfac in cmfacs:
    infilename  = folder+'basHH_cmf'+str(cmfac)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+changestring+'_V_sprxf.txt' 
    infile = open(infilename,'r')
    lines = infile.readlines()
    N = len(lines) # Do I need this?
    
    time = []
    V    = []
    Ca   = []
    
    for line in lines:
        words = line.split()
        time.append(float(words[0]))
        V.append(float(words[1]))
    infile.close()
    
    time = np.array(time)
    V    = np.array(V)

    idxs = np.where(time<(idelay+idur-1))
    li   = idxs[-1] # Want the last one
    
    Vmax = V[0] # max(V) might appear at the end, which is of no use to us
    Vmin = min(V[li])
    dV = Vmax-Vmin
    dV1e = dV*math.exp(-1)
    V1e  = Vmin+dV1e
    print('Vmax:',Vmax)
    print('Vmin:',Vmin)
    print('V1e:',V1e)

    # Will perform linear interpolation near the crossing with V1e.
    tbef = 0
    taft = 0
    Vbef = 0
    Vaft = 0
    for i in range(len(V)):
        if V[i]-V1e<0:
            tbef = time[i-1]
            taft = time[i]
            Vbef = V[i-1]
            Vaft = V[i]
            print('V:',V[i-1],'i-1:',i-1,'; t:', time[i-1])
            print('V:',V[i],'i:',i,'; t:', time[i])
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
        plt.plot(time,V, label='Data')
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
        plt.axis([idelay-0.5,idelay+10,Vmin-0.05,Vmax+0.05])
        plt.legend(loc='upper right')
        plt.show()
        
    tau = tau3
    
    # Conversions (to SI):
    iamp_SI=iamp*1e-9    # Current now in Amperes
    tau*=1e-3            # Time now in seconds
    
    dV*= 1e-3            # Voltage now in Volts
    R  = abs(dV/iamp_SI) # Resistance in Ohms
    C  = tau/R           # Capacitance in Farads
    C *= 1e12            # Capacitance in picoFarads
    #Cm = C/A             # Specific capacitance. C in Farads, A in (mu m)^2
    #Cm*= 1e14            # Conversion to get muF/cm^2
    
    tau*=1e3             # Now in ms again.
    
    outcms.append(C)
    outtaus.append(tau)
    Rins.append(R)
    outfile.write('%.2f %.12f\n' % (cmfac,C))
    outfile_R.write('%.2f %.12f\n' % (cmfac,R))
    outfile_tau.write('%.2f %.12f\n' % (cmfac,tau))
outfile.close()
outfile_R.close()
outfile_tau.close()

plt.figure(figsize=(6,5))
plt.plot(cmfacs,outcms,'-o')
plt.xlabel(r'$C_m$ [$\mu$F/cm$^2$] (Cell parameter)')
plt.ylabel(r'$C$ [$pF$] (Measured value)')
plt.title(r'Measured $C_m$ vs cell parameter $C_m$')
plt.tight_layout()
plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.plot(cmfacs,outtaus,'-o')
plt.xlabel(r'$C_m$ [$\mu$F/cm$^2$] (Cell parameter)')
plt.ylabel(r'$\tau_m$ [ms] (Measured value)')
plt.title(r'Measured $\tau_m$ vs cell parameter $C_m$')
plt.tight_layout()
plt.savefig(plotname_tau)

plt.figure(figsize=(6,5))
plt.plot(cmfacs,Rins,'-o')
plt.xlabel(r'$C_m$ [$\mu$F/cm$^2$] (Cell parameter)')
plt.ylabel(r'$R_{in}$ [$\Omega$] (Measured value)')
plt.title(r'Input resistance $R$ vs cell parameter $C_m$')
plt.tight_layout()
plt.savefig(plotname_R)

print('varyE:',varyE)
print('varyg:',varyg)