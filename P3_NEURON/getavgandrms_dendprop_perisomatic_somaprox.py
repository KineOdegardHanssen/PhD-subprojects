import matplotlib.pyplot as plt
import numpy as np

def avg_and_rms(x):
    N = len(x)
    avgx = np.mean(x)
    print('avgx:',avgx)
    rmsx = 0
    for i in range(N):
        rmsx += (avgx-x[i])**2
    rmsx = np.sqrt(rmsx/(N-1))
    return avgx,rmsx

plottitletext = 'soma and prox. dend.'

##### Adjustable parameters/specifications #########
iduraa = 2
# File/simulation selection:
testmodels = [478513437,478513407,488462965] 
Nmodels = len(testmodels)
idur = 2 # ms
iamp = 1.0 # nA

varycm = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0]
N = len(varycm)
propvels = np.zeros(N)
propvel_rmss = np.zeros(N)

outfolder = 'Cell_average/'
cellmodelstring = ''
for i in range(Nmodels):
    cellmodelstring = cellmodelstring +str(testmodels[i])+'_' # Will get a trailing underscore, but whatever
outfolder = outfolder + cellmodelstring +'/'

outfilename_avg = outfolder+'idur%i_iamp' % idur+str(iamp)+ '_somaprox_addedRa_avgrmsallpropvel_ends.txt'
plotname = outfolder +'idur%i_iamp' % idur+str(iamp)+ '_somaprox_addedRa_avgrmsallpropvel_ends.png'
outfile = open(outfilename_avg,'w')

for i in range(N):
    propvels_all = []
    for k in range(Nmodels):
        testmodel = testmodels[k]
        if testmodel==496497595:
            cm_soma = 1.14805
            cm_dend = 9.98231
            cm_axon = 3.00603
            v_init = -86.5 # mV
        elif testmodel==488462965:
            cm_soma = 3.31732779736
            cm_dend = 3.31732779736
            cm_axon = 3.31732779736
            v_init = -86.5 # mV
        elif testmodel==480633479:
            cm_soma = 0.704866 # 0.704866118957
            cm_dend = 0.704866 # 0.704866118957
            cm_axon = 0.704866 # 0.704866118957
            v_init = -96.8
        elif testmodel==478513407:
            cm_soma = 1.0
            cm_dend = 1.0
            cm_axon = 1.0
            v_init = -83.7
        elif testmodel==478513437:
            cm_soma = 2.34539964752
            cm_dend = 2.34539964752
            cm_axon = 2.34539964752
            v_init = -86.8
        # Update cm
        cm_soma = varycm[i]
        
        folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/dendritepropagation/'
        # Update file name
        infilename = folder+'idur%i_iamp' % idur+str(iamp)+'_cmspr'+str(cm_soma) + '_cmd'+str(cm_dend)+'_cma'+str(cm_axon)+'_vi'+str(v_init)+ '_wRa_somaprox_allpropvels_ends.txt' # Needed to cut down on the name length
        print('infilename:',infilename)
        infile = open(infilename,'r')
        lines = infile.readlines()
        for line in lines:
            propvel = float(line.split()[0])
            propvels_all.append(propvel)
        infile.close()
    propvels[i],propvel_rmss[i]=avg_and_rms(propvels_all)
    outfile.write('%.2f %.16f %.16f\n' % (varycm[i],propvels[i],propvel_rmss[i]))
outfile.close()

# Plot results
plt.figure(figsize=(6,5))
plt.errorbar(varycm,propvels,yerr=propvel_rmss,capsize=2)
plt.xlabel(r'$C_{m}$ of %s [$\mu$ F/cm$^2$]' % plottitletext)
plt.ylabel(r'Signal velocity along dendrites [m/s]')
plt.title(r'Signal velocity vs %s $C_{m}$' % plottitletext)
plt.tight_layout()
plt.savefig(plotname)
