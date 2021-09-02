import numpy as np
import matplotlib.pyplot as plt
import math

dtexp = -8

varydiam = False # True # 

zoomed = False
somasize = 10
dendlen  = 1000 # 1 # 2 # 5 # 10 # 20 # 50 # 100 # 200 # 500 #
denddiam = 2 # 0.01 # 20 # 5 # 10 # 20 # 
if varydiam==True:
    denddiams = [0,0.01,0.1,1,2,5,10,20]
    N = len(denddiams)
else:
    dendlens = [0,1,2,5,10,20,50,100,200,500,1000]
    N = len(dendlens)
    
A      = np.pi*somasize*somasize
Ra     = 100
v_init = -65
gpas   = 0.0003
vpas   = -65
epas   = vpas

idur = 100   # ms
iamp = -0.1  # nA 
currentfolder = 'current_idur%i_iamp'%idur+str(iamp)+'/dtexp%i/' % dtexp

taus_at_1p5  = []
filenames    = []
Rins_at_1p5  = []
Cs_at_1p5    = []
filenames_Cm = []
filenames_R  = []
measured_C   = []
measured_R   = []
somaonly_folder = 'Somaonly/pas/Results/IStim/Soma10/'+currentfolder
filename_nodend = somaonly_folder+'somaonly_cms_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_gpas'+str(gpas)+'_epas'+str(epas)+'_tau_highCms.txt'
filename_nodend_Cm  = somaonly_folder+'somaonly_cms_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_gpas'+str(gpas)+'_epas'+str(epas)+'_Cm_highCms.txt'
filename_nodend_R = somaonly_folder+'somaonly_cms_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_gpas'+str(gpas)+'_epas'+str(epas)+'_Rin_highCms.txt'
filenames.append(filename_nodend)
filenames_Cm.append(filename_nodend_Cm)
filenames_R.append(filename_nodend_R)
plotfolder = 'Comparemodels/BAS_vs_somaonly_passive/dtexp%i/' % dtexp

if varydiam==True:
    endsnippet = '_dtexp%i_varydenddiam' % dtexp
    otherparamsn = '_len%i' % dendlen
else:
    endsnippet = '_dtexp%i_varydendlen' % dtexp
    otherparamsn = '_diam' + str(denddiam)

plotname   = plotfolder + 'BAS_vs_onecomp_cms_idur%i_iamp'%idur+str(iamp)+'_Ra%i_gpas'%Ra+str(gpas)+'_vpas' +str(vpas)+otherparamsn+'+endsnippet
plotname_2 = plotfolder + 'BAS_vs_onecomp_cm1.5_idur%i_iamp'%idur+str(iamp)+'_Ra%i_gpas'%Ra+str(gpas)+'_vpas' +str(vpas)+otherparamsn+'_tau'+endsnippet+'_percentagediff.png'
plotname_C = plotfolder + 'BAS_vs_onecomp_cm1.5_idur%i_iamp'%idur+str(iamp)+'_Ra%i_gpas'%Ra+str(gpas)+'_vpas' +str(vpas)+otherparamsn+'_tau'+endsnippet+'_measureC.png'
plotname_R = plotfolder + 'BAS_vs_onecomp_cm1.5_idur%i_iamp'%idur+str(iamp)+'_Ra%i_gpas'%Ra+str(gpas)+'_vpas' +str(vpas)+otherparamsn+'_tau'+endsnippet+'_Rin.png'
plotname_tau_vs_d = plotfolder + 'BAS_vs_onecomp_cm1_idur%i_iamp'%idur+str(iamp)+'_Ra%i_gpas'%Ra+str(gpas)+'_vpas' +str(vpas)+otherparamsn+'_tau_vs_d'+endsnippet+'.png'
plotname_R_vs_d = plotfolder + 'BAS_vs_onecomp_cm1.5_idur%i_iamp'%idur+str(iamp)+'_Ra%i_gpas'%Ra+str(gpas)+'_vpas' +str(vpas)+otherparamsn+'_tau'+endsnippet+'_Rin_vs_d.png'
plotname_C_vs_d = plotfolder + 'BAS_vs_onecomp_cm1.5_idur%i_iamp'%idur+str(iamp)+'_Ra%i_gpas'%Ra+str(gpas)+'_vpas' +str(vpas)+otherparamsn+'_tau'+endsnippet+'_C_vs_d.png'
if zoomed==True:
    plotname=plotname+'_zoomed.png'
else:
    plotname=plotname+'.png'

for i in range(1,N):
    if varydiam==True:
        denddiam = denddiams[i]
    else:
        dendlen = dendlens[i]
    folder = 'Ball-and-stick models/BAS_passive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+currentfolder
    filename = folder + 'baspass_cms_idur%i_iamp' % idur+str(iamp)+'_Ra%i_gpas'%Ra+str(gpas)+'_vpas' +str(vpas)+'_tau_highCms.txt'
    filename_Cm = folder +'baspass_cms_idur%i_iamp' % idur+str(iamp)+'_Ra%i_gpas'%Ra+str(gpas)+'_vpas' +str(vpas)+'_Cm_highCms.txt'
    filename_R = folder +'baspass_cms_idur%i_iamp' % idur+str(iamp)+'_Ra%i_gpas'%Ra+str(gpas)+'_vpas' +str(vpas)+'_Rin'
    filenames.append(filename)
    filenames_R.append(filename_R)
    filenames_Cm.append(filename_Cm)

plt.figure(figsize=(6,5))
for i in range(N):
    file = open(filenames[i],'r')
    lines = file.readlines()
    
    cms = []
    taus = []
    #print('NUMBER OF LINES:', len(lines), filenames[i]) ### FOR DEBUGGING
    for line in lines:
        words = line.split()
        if len(words)>1:
            cm = float(words[0])
            tau = float(words[1])
            cms.append(cm)
            taus.append(tau)
            if cm==1.5:
                taus_at_1p5.append(tau)
    if zoomed==False:
        if varydiam==True:
            plt.plot(cms,taus,label='Dendrite diam: %.2f' % denddiams[i])
        else:
            plt.plot(cms,taus,label='Dendrite len: %i' % dendlens[i])
    else:
        if varydiam==True:
            plt.plot(cms[2:7],taus[2:7],label='Dendrite diam: %.2f' % denddiams[i])
        else:
            plt.plot(cms[2:7],taus[2:7],label='Dendrite len: %i' % dendlens[i])
    file.close()
    
    filename_Cm = filenames_Cm[i]
    file_Cm = open(filename_Cm,'r')
    lines   = file_Cm.readlines()
    cms     = []
    theseC  = []
    for line in lines:
        words = line.split()
        if len(words)>0:
            cm = float(words[0])
            C = float(words[1])*A
            cms.append(cm)
            theseC.append(C)
            if cm==1.5:
                Cs_at_1p5.append(C)
    measured_C.append(theseC)
    file_Cm.close()
    
    file_R = open(filenames_R[i],'r')
    lines = file_R.readlines()
    cms     = []
    theseR  = []
    for line in lines:
        words = line.split()
        if len(words)>1:
            cm = float(words[0])
            R = float(words[1])
            cms.append(cm)
            theseR.append(R)
            if cm==1.5:
                Rins_at_1p5.append(R)
    measured_R.append(theseR)
    file_R.close()

'''
print('measured_C:',measured_C)
print('measured_R:',measured_R)
print('measured_C[0]:',measured_C[0])
print('measured_R[0]:',measured_R[0])
'''

plt.xlabel(r'$C_m$ [$\mu$F/cm$^2$] (Cell parameter)')
plt.ylabel(r'$\tau_m$ [ms] (Measured value)')
plt.title(r'$\tau_m$ vs $C_m$, passive cell, ball-and-stick')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig(plotname)

fracs = np.zeros(N-1)
fracs12 = np.zeros(N)
fracs15 = np.zeros(N)
plt.figure(figsize=(6,5))
plt.ylabel('% difference')
if varydiam==True:
    ############################## Plotting vs diam ###################################################
    print('dendlen=', dendlen, ', varydiam:')
    print('Difference, one compartment and thickest dendrite (tau0/tau_thick):', taus_at_1p5[0]/taus_at_1p5[-1])
    print('Relative difference:', (taus_at_1p5[0]-taus_at_1p5[-1])/taus_at_1p5[0])
    print('thickness 1:', denddiams[0], 'thickness 2:', denddiams[-1])
    for i in range(1,N):
        #print('i:',i,'; N:',N)
        #print('len(taus_at_1p5):',len(taus_at_1p5))
        #print('len(fracs):',len(fracs))
        fracs[i-1] = 100*(taus_at_1p5[0]-taus_at_1p5[i])/taus_at_1p5[0] # Will be 0 at i=0
    plt.plot(denddiams[1:],fracs)
    plt.xlabel(r'Dendrite diameter ($\mu$m)')
    plt.title(r'Difference, $\tau$, one comp. and BAS, dend.len. %i' % dendlen)
    plt.tight_layout()
    plt.savefig(plotname_2)
    # 
    plt.figure(figsize=(6,5))
    for i in range(N):
        plt.plot(cms,measured_C[i],label=r'denddiam=%s' % str(denddiams[i]))
    plt.xlabel(r'Cell parameter $C_m$')
    plt.ylabel(r'Measured $C$')
    plt.title(r'Measured $C$ vs parameter $C_m$, dend.len. %.2f' % dendlen)
    plt.tight_layout()
    plt.legend()
    plt.savefig(plotname_C)
    #
    plt.figure(figsize=(6,5))
    for i in range(N):
        plt.plot(cms,measured_R[i],label=r'denddiam=%s' % str(denddiams[i]))
    plt.xlabel(r'Cell parameter $C_m$')
    plt.ylabel(r'$R_{in}$ [$\Omega$]')
    plt.title(r'$R_{in}$ vs parameter $C_m$, dend.len. %.2f' % dendlen)
    plt.tight_layout()
    plt.legend()
    plt.savefig(plotname_R)
    # 
    plt.figure(figsize=(6,5))
    plt.plot(denddiams,taus_at_1p5)
    plt.xlabel(r'Dendrite diameter ($\mu$m)')
    plt.ylabel(r'$\tau$ [ms]')
    plt.title(r'$\tau$ vs $d$, dend.len. %i, $C_m=1.5$' % dendlen)
    plt.tight_layout()
    plt.savefig(plotname_tau_vs_d)
    # 
    plt.figure(figsize=(6,5))
    plt.plot(denddiams,Rins_at_1p5)
    plt.xlabel(r'Dendrite diameter ($\mu$m)')
    plt.ylabel(r'$R_{in}$ [$\Omega$]')
    plt.title(r'$R_{in}$ vs $d$, dend.len. %i, $C_m=1.5$' % dendlen)
    plt.tight_layout()
    plt.savefig(plotname_R_vs_d)
    # 
    plt.figure(figsize=(6,5))
    plt.plot(denddiams,Cs_at_1p5)
    plt.xlabel(r'Dendrite diameter ($\mu$m)')
    plt.ylabel('Capacitance')
    plt.title(r'$C$ vs $d$, dend.len. %i, $C_m=1.5$' % dendlen)
    plt.tight_layout()
    plt.savefig(plotname_C_vs_d)
else:
    ############################## Plotting vs len ####################################################
    print('denddiam=', denddiam, ', varylen:')
    print('Difference, one compartment and longest dendrite (tau0/tau_long):', taus_at_1p5[0]/taus_at_1p5[-1])
    print('Relative difference:', (taus_at_1p5[0]-taus_at_1p5[-1])/taus_at_1p5[0])
    print('length 1:', dendlens[0], 'length 2:', dendlens[-1])
    for i in range(1,N):
        #print('i:',i,'; N:',N)
        #print('len(taus_at_1p5):',len(taus_at_1p5))
        #print('len(fracs):',len(fracs))
        fracs[i-1] = 100*(taus_at_1p5[0]-taus_at_1p5[i])/taus_at_1p5[0] # Will be 0 at i=0
    plt.plot(dendlens[1:],fracs)
    plt.xlabel(r'Dendrite length ($\mu$m)')
    plt.title(r'Difference, $\tau$, one comp. and BAS, dend.diam. %.2f' % denddiam)
    plt.tight_layout()
    plt.savefig(plotname_2)
    # 
    plt.figure(figsize=(6,5))
    for i in range(N):
        plt.plot(cms,measured_C[i],label=r'dendlen=%s' % str(dendlens[i]))
    plt.xlabel(r'Cell parameter $C_m$')
    plt.ylabel(r'Measured $C$')
    plt.title(r'Measured $C$ vs parameter $C_m$, dend.diam. %.2f' % denddiam)
    plt.tight_layout()
    plt.legend()
    plt.savefig(plotname_C)
    #
    plt.figure(figsize=(6,5))
    for i in range(N):
        plt.plot(cms,measured_R[i],label=r'dendlen=%s' % str(dendlens[i]))
    plt.xlabel(r'Cell parameter $C_m$')
    plt.ylabel(r'$R_{in}$ [$\Omega$]')
    plt.title(r'$R_{in}$ vs parameter $C_m$, dend.diam. %.2f' % denddiam)
    plt.tight_layout()
    plt.legend()
    plt.savefig(plotname_R)
    #
    plt.figure(figsize=(6,5))
    plt.plot(dendlens,taus_at_1p5)
    plt.xlabel(r'Dendrite length ($\mu$m)')
    plt.ylabel(r'$\tau$ [ms]')
    plt.title(r'$\tau$ vs $l$, dend.diam. %.2f, $C_m=1.5$' % denddiam)
    plt.tight_layout()
    plt.savefig(plotname_tau_vs_d)
    #
    plt.figure(figsize=(6,5))
    plt.plot(dendlens,Rins_at_1p5)
    plt.xlabel(r'Dendrite length ($\mu$m)')
    plt.ylabel(r'$R_{in}$ [$\Omega$]')
    plt.title(r'$R_{in}$ vs $l$, dend.diam. %.2f, $C_m=1.5$' % denddiam)
    plt.tight_layout()
    plt.savefig(plotname_R_vs_d)
    #
    plt.figure(figsize=(6,5))
    plt.plot(dendlens,Cs_at_1p5)
    plt.xlabel(r'Dendrite length ($\mu$m)')
    plt.ylabel('Capacitance')
    plt.title(r'$C$ vs $l$, dend.diam. %.2f, $C_m=1.5$' % (denddiam))
    plt.tight_layout()
    plt.savefig(plotname_C_vs_d)
plt.show()
