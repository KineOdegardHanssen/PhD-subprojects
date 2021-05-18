import matplotlib.pyplot as plt
import numpy as np

##### Adjustable parameters/specifications #########
smallCm = True
iduraa = 2
# File/simulation selection:
testmodel = 488462965#478513407#478513437# 
idur = 2 # ms
iamp = 1.0 # nA
v_init = -86.5 # mV
dendnr = 6

if testmodel==496497595:
    cm_soma = 1.14805
    cm_dend = 9.98231
    cm_axon = 3.00603
elif testmodel==488462965:
    cm_soma = 3.31732779736
    cm_dend = 3.31732779736
    cm_axon = 3.31732779736
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

varycm = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0]
N = len(varycm)

plottitletext = 'soma and proximal dendrites'


folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/dendritepropagation/'
outfolder = folder
outfilename_avg = outfolder+'idur%i_iamp' % idur+str(iamp)+'_cmspr' + '_'+str(testmodel)+ '_'+ '_vinit'+str(v_init)+'_addedRa_somaprox_avgpropvel_ends.txt'
plotname = outfolder +'idur%i_iamp' % idur+str(iamp)+'_cmspr' + '_'+str(testmodel)+ '_'+'_vinit'+str(v_init)+'_addedRa_somaprox_avgpropvel_ends.png'

propvels = np.zeros(N)
propvel_rmss = np.zeros(N)

outfile = open(outfilename_avg,'w')

for i in range(N):
    # Update cm
    cm_soma = varycm[i]
    # Update file name
    if testmodel==496497595:
        infilename = folder+'idur%i_iamp' % idur+str(iamp)+'_cmspr'+str(cm_soma) + '_cmd'+str(cm_dend)+'_cma'+str(cm_axon)+'_vi'+str(v_init)+ '_wRa_somaprox_avgpropvel_ends.txt' # Needed to cut down on the name length
    else:
        infilename = folder+'idur%i_iamp' % idur+str(iamp)+'_cmspr'+str(cm_soma) + '_cmd'+str(cm_dend)+'_cma'+str(cm_axon)+'_vi'+str(v_init)+ '_wRa_somaprox_avgpropvel_ends.txt' # Needed to cut down on the name length
    print('infilename:',infilename)
    print('outfilename:',outfilename_avg)
    infile = open(infilename,'r')
    line = infile.readline()
    print('i:',i, '; cm_soma:',cm_soma)
    print('Line:', line)
    propvel = float(line.split()[0])
    propvel_rms = float(line.split()[1])
    propvels[i] = propvel
    propvel_rmss[i] = propvel_rms
    outfile.write('%.11f %.16f %.16f\n' % (varycm[i],propvel,propvel_rms))
    infile.close()
outfile.close()

# Plot results
plt.figure(figsize=(6,5))
plt.errorbar(varycm,propvels,yerr=propvel_rmss,capsize=2)
plt.xlabel(r'$C_{m}$ of %s [$\mu$ F/cm$^2$]' % plottitletext)
plt.ylabel(r'Signal velocity along dendrites [m/s]')
plt.title(r'Signal velocity vs %s $C_{m}$' % plottitletext)
plt.tight_layout()
plt.savefig(plotname)