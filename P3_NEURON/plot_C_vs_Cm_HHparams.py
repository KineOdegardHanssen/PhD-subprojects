import numpy as np
import matplotlib.pyplot as plt
import math

dtexp = -8


zoomed = False
somasize = 10
dendlen  = 1000 
denddiam = 2 
varymech = 'Na' # 'K' # 'leak'
varyE_bool = True
varyEs = [30,40,50,60,70]
varygs = 'None'
    
varylist = [] # Should be redundant
plotstring = '_vary'
if varyE_bool==True:
    #varylist = varyE
    plotstring = plotstring + 'E'
    N = len(varyEs)
else:
    #varylist = varyg
    plotstring = plotstring + 'g'
    N = len(varygs)
Nvary    = len(varylist)

if varymech=='Na':
    varyfolder = 'VaryNa/'
elif varymech=='K':
    varyfolder = 'VaryK/'
else:
    varyfolder = 'VaryLeak/'

changestring =''
if varyE_bool==True:
    changestring = changestring+'_varyE_gdflt'
else:
    changestring = changestring+'_Edefault_varyg'

A      = np.pi*somasize*somasize
Ra     = 100
v_init = -65
gpas   = 0.0003
vpas   = -65
epas   = vpas

idur = 1000   # ms
iamp = -0.1  # nA 
currentfolder = 'current_idur%i_iamp'%idur+str(iamp)+'/'

filenames    = []
filenames_Cm = []
filenames_R  = []
measured_C   = []
measured_R   = []
plotfolder = 'Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+varyfolder+currentfolder


plotname_C   = plotfolder + 'BASHH_cmfs_idur%i_iamp'%idur+str(iamp)+'_Ra%i' %Ra+changestring +'_CvsCm'
plotname_R   = plotfolder + 'BASHH_cmfs_idur%i_iamp'%idur+str(iamp)+'_Ra%i' %Ra+changestring +'_RinvsCm'
plotname_tau = plotfolder + 'BASHH_cmfs_idur%i_iamp'%idur+str(iamp)+'_Ra%i' %Ra+changestring +'_tauvsCm'
if zoomed==True:
    plotname_C   = plotname_C+'_zoomed.png'
    plotname_R   = plotname_R+'_zoomed.png'
    plotname_tau = plotname_tau+'_zoomed.png'
else:
    plotname_C = plotname_C+'.png'
    plotname_R = plotname_R+'.png'
    plotname_tau = plotname_tau+'.png'

for i in range(N):
    changestring =''
    if varyE_bool==True:
        varyE = varyEs[i]
        changestring = changestring+'_E'+str(varyE)+'_gdflt'
    else:
        varyg = varygs[i]
        changestring = changestring+'_Edefault_g'+str(varyg)

    folder = 'Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+varyfolder+currentfolder
    filename    = folder + 'basHH_cmfs_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+changestring+'_tau_highCms.txt'
    filename_Cm = folder + 'basHH_cmfs_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+changestring+'_Ctot_highCms.txt'
    filename_R  = folder +'basHH_cmfs_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+changestring+'_Rin'#_highCms.txt'
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
    if zoomed==False:
        if varyE_bool==True:
            plt.plot(cms,taus,label='E=%i' % varyEs[i])
        else:
            plt.plot(cms,taus,label='g=%i' % varygs[i])
    else:
        if varyE_bool==True:
            plt.plot(cms[2:7],taus[2:7],label='E=%i' % varyEs[i])
        else:
            plt.plot(cms[2:7],taus[2:7],label='g=%i' % varygs[i])
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
            C = float(words[1])
            cms.append(cm)
            theseC.append(C)
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
plt.title(r'$\tau_m$ vs $C_m$, ball-and-stick HH, l=%i, d=%.2f, HH vary %s' % (dendlen,denddiam,varymech))
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig(plotname_tau)

if varyE_bool==True:
    ############################## Plotting vs E ###################################################
    # 
    plt.figure(figsize=(6,5))
    for i in range(N):
        plt.plot(cms,measured_C[i],label=r'E=%s' % str(varyEs[i]))
    plt.xlabel(r'Cell parameter $C_m$ ($\mu$F/cm$^2$)')
    plt.ylabel(r'Measured $C$ (pF)')
    plt.title(r'$C$ vs parameter $C_m$, l=%i, d=%.2f, HH vary %s' % (dendlen,denddiam,varymech))
    plt.tight_layout()
    plt.legend()
    plt.savefig(plotname_C)
    #
    plt.figure(figsize=(6,5))
    for i in range(N):
        plt.plot(cms,measured_R[i],label=r'E=%s' % str(varyEs[i]))
    plt.xlabel(r'Cell parameter $C_m$ ($\mu$F/cm$^2$)')
    plt.ylabel(r'$R_{in}$ [$\Omega$]')
    plt.title(r'$R_{in}$ vs parameter $C_m$, l=%i, d=%.2f, HH vary %s' % (dendlen,denddiam,varymech))
    plt.tight_layout()
    plt.legend()
    plt.savefig(plotname_R)
else:
    ############################## Plotting vs g ####################################################    
    plt.figure(figsize=(6,5))
    for i in range(N):
        plt.plot(cms,measured_C[i],label=r'g=%s' % str(varygs[i]))
    plt.xlabel(r'Cell parameter $C_m$ ($\mu$F/cm$^2$)')
    plt.ylabel(r'Measured $C$ (pF)')
    plt.title(r'$C$ vs parameter $C_m$, l=%i, d=%.2f, HH vary %s' % (dendlen,denddiam,varymech))
    plt.tight_layout()
    plt.legend()
    plt.savefig(plotname_C)
    #
    plt.figure(figsize=(6,5))
    for i in range(N):
        plt.plot(cms,measured_R[i],label=r'g=%s' % str(varygs[i]))
    plt.xlabel(r'Cell parameter $C_m$ ($\mu$F/cm$^2$)')
    plt.ylabel(r'$R_{in}$ [$\Omega$]')
    plt.title(r'$R_{in}$ vs parameter $C_m$, l=%i, d=%.2f, HH vary %s' % (dendlen,denddiam,varymech))
    plt.tight_layout()
    plt.legend()
    plt.savefig(plotname_R)
    #
plt.show()
