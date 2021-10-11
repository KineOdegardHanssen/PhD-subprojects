import numpy as np
import matplotlib.pyplot as plt
import math

dtexp = -8

testit = False # True #  #If I am unsure about the results, I can check the fit.

varymech = 'Na' # 'K' # 'leak'
varyE_bool = True
varyE = 50 #[-10,20,50,60]
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

changestring =''
if varyE_bool==True:
    varyE = elem
    changestring = changestring+'_E'+str(varyE)+'_gdflt'
else:
    varyg = elem
    changestring = changestring+'_Edefault_g'+str(varyg)

if varymech=='NaV':
    folderstring = 'VaryNa/' 
    plotstring = plotstring + '_NaV'
elif varymech=='pas':
    folderstring = 'VaryPas/'
    plotstring = plotstring + '_Pas'
elif varymech=='Kd':
    folderstring = 'VaryKd/'
    plotstring = plotstring + '_Kd'
elif varymech=='Kv2like':
    folderstring = 'VaryKv2like/'
    plotstring = plotstring + '_Kv2like'
elif varymech=='Kv3_1':
    folderstring = 'VaryKv3_1/'
    plotstring = plotstring + '_Kv3_1'
elif varymech=='SK':
    folderstring = 'VarySK/'
    plotstring = plotstring + '_SK'
elif varymech=='K_T':
    folderstring = 'VaryK_T/'
    plotstring = plotstring + '_K_T'
elif varymech=='Im_v2':
    folderstring = 'VaryIm_v2/'
    plotstring = plotstring + '_Im_v2'

# loop de loop? # I really should...
somasize = 10
dendlen  = 1000
denddiam = 0.01
L      = somasize
diam   = somasize
A      = np.pi*L*diam
Ra     = 100
v_init = -65
gpas   = 0.0003
vpas   = -65
long   = True

models = [485694403,478513407,478513437,488462965,489931686] # Order...
Nmodels = len(models)
    
# Values of membrane capacitance (parameters):
cmfacs = [0.1,0.5,0.75,0.88,0.94,1.0,1.25,1.5,2.0,3.0]
incms_all   = []
outcs_all   = []
outtaus_all = [] # I'm not using these yet
Rins_all    = [] # I'm not using these yet

# Change current:
idur = 100   # ms
iamp = -0.1  # nA 
idelay = 10

if long==True:
    plotsuffix = '_highCms.png'
    filesuffix = '_highCms.txt'
else:
    plotsuffix = '.png'
    filesuffix = '.txt'

plotfolder   = 'Results/Soma%i/' % somasize
plotname_all = plotfolder +  'baspv_cmfs_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+changestring+'_sprx+'_Ctot_all'+plotsuffix
plotname_avg = plotfolder +  'baspv_cmfs_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+changestring+'_sprx+'_Ctot_avgrms'+plotsuffix
for model in models:
    folder = 'Results/%i/IStim/Soma%i/dendlen%i/denddiam'% (model,somasize,dendlen)+str(denddiam) +'/current_idur%i_iamp'% idur+str(iamp)+'/'
    infilename     = folder +'baspv_cmfs_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+changestring+'_sprx+'_Ctot'+filesuffix
    
    incms = []
    outcs = []
    
    infile = open(infilename,'r')
    lines  = infile.readlines()
    
    for line in lines:
        words = line.split()
        if len(words)>0:
            incms.append(float(words[0]))
            outcs.append(float(words[1]))
    infile.close()
    incms_all.append(incms)
    outcs_all.append(outcs)
    Ncms = len(incms) # Assuming that we have the same number of cms for all models

outcs_avg = np.zeros(Ncms)
outcs_rms = np.zeros(Ncms)

for i in range(Nmodels):
    outcs = outcs_all[i]
    for j in range(Ncms):
        outcs_avg[j]+=outcs[j]
outcs_avg/=Nmodels

for i in range(Nmodels):
    outcs = outcs_all[i]
    for j in range(Ncms):
        outcs_rms[j]+=(outcs[j]-outcs_avg[j])**2
for i in range(Ncms):
    outcs_rms[i]=np.sqrt(outcs_rms/(Nmodels-1))

plt.figure(figsize=(6,5))
plt.errorbar(incms,outcs_avg,yerr=outcs_rms,capsize=2)
plt.xlabel(r'$C_m$ [$\mu$F/cm$^2$] (Cell parameter)')
plt.ylabel(r'$C$ [$pF$] (Measured value)')
plt.title(r'Measured $C$ vs cell parameter $C_m$')
plt.savefig(plotname_avgrms)

plt.figure(figsize=(6,5))
for i in range(Nmodels):
    plt.plot(incms,outcs[i],label='%i' % models[i])
plt.xlabel(r'$C_m$ [$\mu$F/cm$^2$] (Cell parameter)')
plt.ylabel(r'$C$ [$pF$] (Measured value)')
plt.title(r'Measured $C$ vs cell parameter $C_m$')
plt.legend(loc='upper left')
plt.savefig(plotname_all)
