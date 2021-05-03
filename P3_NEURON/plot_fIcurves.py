import numpy as np
import matplotlib.pyplot as plt

testmodels = [478513437,478513407,488462965]
idur   = 1000 # ms
idelay = 100
v_init = -86.5 # mV
idur_SI = idur/1000. # Time is measured in ms
Nmodels = len(testmodels)

iamps = [0,0.02,0.04,0.06,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,
0.38,0.4,0.42,0.44,0.46,0.48,0.5]
Ni    = len(iamps)
fs_store     = []
fs_rms_store = []
    
# Defaulting to original values:
# DO NOT TOUCH THESE!
# SET THEM BELOW INSTEAD!
cm_somas = [0.5,1.0,2.0] #[0.5,1.0,1.25,2.0]
Ncm = len(cm_somas)


outfolder = 'Comparemodels/'
cellmodelstring = ''
for i in range(Nmodels):
    cellmodelstring = cellmodelstring +str(testmodels[i])+'_' # Will get a trailing underscore, but whatever
outfolder = outfolder + cellmodelstring +'/'
plotname  = outfolder+'f_I_curves_avg_and_rms.png' 

folder = 'Allen_test_changecapacitance/figures/'
for i in range(Ncm):
    cm_soma = cm_somas[i]
    fs     = np.zeros(Ni)
    fs_rms = np.zeros(Ni)
    f_thisCm = np.zeros((Nmodels,Ni))
    k = 0
    for testmodel in testmodels:
        infolder   = folder+'%i/' % testmodel
        infilename = infolder+'f_I_curves/f_I_curves_cmsoma'+str(cm_soma)+'.txt'
        infile = open(infilename,'r')
        lines = infile.readlines()
        
        j = 0
        for line in lines:
            words = line.split()
            iamp_this = float(words[0])
            f_this    = float(words[1])
            print('iamp:',iamp_this,'; f:',f_this,'; model:',testmodel)
            fs[j] += f_this
            f_thisCm[k,j]=f_this
            print('fs[',j,']:',fs[j])
            j+=1
        k+=1
    fs/=Nmodels
    print('fs:',fs)
    for l in range(Ni):
        for m in range(Nmodels):
            fs_rms[l]+=(fs[l]-f_thisCm[m,l])**2
        fs_rms[l]=np.sqrt(fs_rms[l]/(Nmodels-1))
    fs_store.append(fs)
    fs_rms_store.append(fs_rms)

plt.figure(figsize=(6,5))
for i in range(Ncm):
    fs = fs_store[i]
    fs_rms = fs_rms_store[i]
    plt.errorbar(iamps[:-5],fs[:-5],yerr=fs_rms[:-5],capsize=2,label=r'$C_m$=%.2f $\mu$F/cm$^2$' % cm_somas[i])
plt.xlabel(r'Input current [nA]')
plt.ylabel(r'Frequency [Hz]')
plt.title(r'f-I curve')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig(plotname)
plt.show()