import numpy as np
import matplotlib.pyplot as plt
import math

dtexp = -8

# loop de loop? # Yes.
L      = 10
diam   = 10
A      = np.pi*L*diam
Ra     = 100
v_init = -65
long   = True # False # 
    
# Values of membrane capacitance (parameters):
cms = [0.1,0.5,0.75,0.88,0.94,1.0,1.25,1.5,2.0,3.0]
outcms  = []
outtaus = []
Rins    = []

# Change current:
idur   = 100   # ms
iamp   = -0.1  # nA 
idelay = 10 
gpas   = 0.0003
epas   = -65
testit = False # True # If I am unsure about the results, I can check the fit.

folder = 'Results/IStim/Soma%i/current_idur%i_iamp'% (L,idur)+str(iamp)+'/dtexp%i/' % dtexp
outfilename = folder +'somaonly_cms_idur%i_iamp' % idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_gpas'+str(gpas)+'_epas'+str(epas)+'_Cm'
outfilename_tau = folder +'somaonly_cms_idur%i_iamp' % idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_gpas'+str(gpas)+'_epas'+str(epas)+'_tau'
outfilename_R = folder +'somaonly_cms_idur%i_iamp' % idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_gpas'+str(gpas)+'_epas'+str(epas)+'_Rin'
if long==True:
    plotname        = outfilename + '_highCms.png'
    plotname_tau    = outfilename_tau + '_highCms.png'
    plotname_R      = outfilename_R + '_highCms.png'
    outfilename     = outfilename + '_highCms.txt'
    outfilename_tau = outfilename_tau + '_highCms.txt'
    outfilename_R   = outfilename_R + '_highCms.txt'
else:
    plotname        = outfilename + '.png'
    plotname_tau    = outfilename_tau + '.png'
    plotname_R      = outfilename_R + '.png'
    outfilename     = outfilename + '.txt'
    outfilename_tau = outfilename_tau + '.txt'
    outfilename_R   = outfilename_R + '.txt'
outfile = open(outfilename,'w')
outfile_R = open(outfilename_R,'w')
outfile_tau = open(outfilename_tau,'w')

for cm in cms:
    infilename  = folder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_gpas'+str(gpas)+'_epas'+str(epas)+'_V.txt' 
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
    
    Vmax = max(V)
    Vmin = min(V)
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
    testit= False # True #
    if testit==True:
        N = 100*idur
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
        plt.axis([idelay-0.5,idelay+10,min(V)-0.01,max(V)+0.01])
        plt.legend(loc='upper right')
        plt.show()
        
    tau = tau3
    
    # Conversions (to SI):
    iamp_SI=iamp*1e-9    # Current now in Amperes
    tau*=1e-3            # Time now in seconds
    
    dV*= 1e-3            # Voltage now in Volts
    R  = abs(dV/iamp_SI) # Resistance in Ohms
    C  = tau/R           # Capacitance in Farads
    Cm = C/A             # Specific capacitance. C in Farads, A in (mu m)^2
    Cm*= 1e14            # Conversion to get muF/cm^2
    
    outcms.append(Cm)
    outtaus.append(tau)
    Rins.append(R)
    outfile.write('%.2f %.12f\n' % (cm,Cm))
    outfile_R.write('%.2f %.12f\n' % (cm,R))
    outfile_tau.write('%.2f %.12f\n' % (cm,tau))
outfile.close()
outfile_R.close()
outfile_tau.close()

print('outfilename_R:',outfilename_R)

plt.figure(figsize=(6,5))
plt.plot(cms,outcms,'-o')
plt.xlabel(r'$C_m$ [$\mu$F/cm$^2$] (Cell parameter)')
plt.ylabel(r'$C_m$ [$\mu$F/cm$^2$] (Measured value)')
plt.title(r'Measured $C_m$ vs cell parameter $C_m$')
plt.tight_layout()
plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.plot(cms,outtaus,'-o')
plt.xlabel(r'$C_m$ [$\mu$F/cm$^2$] (Cell parameter)')
plt.ylabel(r'$\tau_m$ [ms] (Measured value)')
plt.title(r'Measured $\tau_m$ vs cell parameter $C_m$')
plt.tight_layout()
plt.savefig(plotname_tau)

plt.figure(figsize=(6,5))
plt.plot(cms,Rins,'-o')
plt.xlabel(r'$C_m$ [$\mu$F/cm$^2$] (Cell parameter)')
plt.ylabel(r'$R_{in}$ [$\Omega$] (Measured value)')
plt.title(r'Input resistance $R$ vs cell parameter $C_m$')
plt.tight_layout()
plt.savefig(plotname_R)