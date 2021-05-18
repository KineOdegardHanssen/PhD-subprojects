import matplotlib.pyplot as plt
import numpy as np


# Varycm: All, soma or dendrite
#varywhichcm = 'a' # All # Only for perisomatic so far
varywhichcm = 's' # Soma
#varywhichcm = 'd' # Dendrite

##### Adjustable parameters/specifications #########
smallCm = True
iduraa = 2
# File/simulation selection:
testmodel = 478513437#488462965#478513407#  ##480633479#496497595#496497595
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

varycm = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
N = len(varycm)

plottitletext = 'soma and prox. dend.'


folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/dendritepropagation/' # Yeah, this works for all.
outfolder = folder
outfilename_avg = outfolder+'idur%i_iamp' % idur+str(iamp)+'_cmspr' + '_'+str(testmodel)+ '_somaprox_vinit'+str(v_init)+'_addedRa_avgVmax_ends.txt'
plotname = outfolder +'idur%i_iamp' % idur+str(iamp)+ '_'+str(testmodel)+ '_somaprox_vinit'+str(v_init)+'_addedRa_avgVmax_ends.png'

vmax_avgs = np.zeros(N)
vmax_rmss = np.zeros(N)

outfile = open(outfilename_avg,'w')

for i in range(N):
    #print('cm:',varycm[i])
    # Update cm
    if varywhichcm=='a':
        cm_soma = varycm[i]
        cm_dend = varycm[i]
        cm_axon = varycm[i]
    elif varywhichcm=='s':
        cm_soma = varycm[i]
    else:
        cm_dend = varycm[i]
    # Update file name
    infilename = folder+'idur%i_iamp' % idur+str(iamp)+'_cmspr' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_somaprox_avgpropvel_Vmax_ends.txt'
    #print('infilename:',infilename)
    #print('outfilename:',outfilename_avg)
    infile = open(infilename,'r')
    line = infile.readline()
    #print('i:',i, '; cm_soma:',cm_soma)
    #print('Line:', line)
    vmax_avg = float(line.split()[0])
    vmax_rms = float(line.split()[1])
    vmax_avgs[i] = vmax_avg
    vmax_rmss[i] = vmax_rms
    outfile.write('%.11f %.16f %.16f\n' % (varycm[i],vmax_avg,vmax_rms))
    infile.close()
outfile.close()

#print('len(varycm):',len(varycm))
#print('len(vmax_avgs):',len(vmax_avgs))
#print('len(vmax_rmss):',len(vmax_rmss))

# Plot results
plt.figure(figsize=(6,5))
plt.errorbar(varycm,vmax_avgs,yerr=vmax_rmss,capsize=2)
plt.xlabel(r'$C_{m}$ of %s [$\mu$ F/cm$^2$]' % plottitletext)
plt.ylabel(r'Maximum membrane potential $V_{max}$ in dendrite ends [mV]')
plt.title(r'$V_{max}$ vs %s $C_{m}$' % plottitletext)
plt.tight_layout()
plt.savefig(plotname)

print('subfolder:', subfolder)