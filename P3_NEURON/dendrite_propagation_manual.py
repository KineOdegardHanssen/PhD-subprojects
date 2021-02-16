import matplotlib.pyplot as plt
import numpy

##### Adjustable parameters/specifications #########
# File/simulation selection:
testmodel = 488462965 # 496497595 #
iduraa = 2
if testmodel==496497595:
    idur = iduraa # ms
elif testmodel==488462965:
    idur = 2 # ms
iamp = 1.0 # nA
v_init = -86.5 # mV

dendnr = 6
L      = 201.927
if testmodel==496497595:
    idxs = [50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66]
elif testmodel==488462965:
    idxs = [32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42]

Nidx   = len(idxs)
dL     = L/float(Nidx)
####################################################

# Defaulting to original values:
# DO NOT TOUCH THESE!
if testmodel==496497595:
    cm_soma = 1.14805
    cm_dend = 9.98231
    cm_axon = 3.00603
elif testmodel==488462965:
    cm_soma = 3.31732779736
    cm_dend = 3.31732779736
    cm_axon = 3.31732779736

# You can touch these instead:
#cm_soma = 3.0
cm_dend = 5.0 #cm_soma #5.0
#cm_axon = cm_soma #3.00603

peaktimes = numpy.zeros(Nidx)


folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/dendritepropagation/'

for j in range(Nidx):
    #### OLD: #######
    #if testmodel==496497595:
    #    filename = folder+'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) +'_vinit'+str(v_init)+ '_addedRa_V_dendsec%i_dend%i.txt' % (j,dendnr)
    #else: # Shorter name due to too many decimals
    #################
    filename = folder+"idur%i_iamp" % idur + str(iamp)+"_cms" + str(cm_soma) + "_cmd" + str(cm_dend) + "_cma"+ str(cm_axon) + "_vinit"+str(v_init)+"_wRa_V_dsec%i_d%i.txt" % (j,dendnr)
    
    times = []
    V     = []

    # Read from file
    infile = open(filename,'r')
    lines   = infile.readlines()
    for line in lines:
        words = line.split()
        if len(words)!=0: # words[0]:time, words[1]:V, words[2]:[Ca^2+]
            times.append(float(words[0]))
            V.append(float(words[1]))
    infile.close()
    
    peaktime = -100
    # Ugh, find the time when V is max instead...
    #for i in range (1,len(V)-1):
    #    if V[i-1]<V[i] and V[i+1]<V[i] and V:
    #        peaktimes[j] = times[i]
    #        break
    max_value = max(V)
    peaktimes[j] = times[V.index(max_value)]

# File and plot names
if testmodel==496497595:
    # Plot name:
    figname = folder+'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) +'_vinit'+str(v_init)+ '_addedRa_V_dend%i_propagation.png' % dendnr
    figname_vel = folder+'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) +'_vinit'+str(v_init)+ '_addedRa_V_dend%i_propvel.png' % dendnr # Needed to cut down on the name length
    # Outfilename:
    outfilename_vel = folder+'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) +'_vinit'+str(v_init)+ '_addedRa_V_dend%i_propvel.txt' % dendnr # Needed to cut down on the name length
else:
    figname = folder+'idur%i_iamp' % idur+str(iamp)+'_cms' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vinit'+str(v_init)+ '_wRa_V_dend%i_prop.png' % dendnr
    figname_vel = folder+'idur%i_iamp' % idur+str(iamp)+'_cms' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_V_dend%i_propvel.png' % dendnr # Needed to cut down on the name length
    # Outfilename:
    outfilename_vel = folder+'idur%i_iamp' % idur+str(iamp)+'_cms' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_V_dend%i_propvel.txt' % dendnr # Needed to cut down on the name length

# Plotting:

plt.figure(figsize=(6,5))
plt.plot(idxs, peaktimes)
plt.xlabel('Index')
plt.ylabel('Peak time [ms]')
plt.title('Peak time for signal along dendrite %i' % dendnr)
plt.tight_layout()
plt.savefig(figname)

# Calculate speed:
vel = numpy.zeros(Nidx-1)
for i in range (0,Nidx-1):
    vel[i] = dL/(peaktimes[i+1]-peaktimes[i])

plt.figure(figsize=(6,5))
plt.plot(idxs[0:Nidx-1], vel)
plt.xlabel('Index')
plt.ylabel('Propagation velocity [$10^{-3}$m/s]')
plt.title('Propagation velocity for signal in each point along dendrite %i' % dendnr)
plt.tight_layout()
plt.savefig(figname_vel)

# Double check:
print('Peaktimes:', peaktimes)

# Other approach: print L/(peaktimes[Nidx-1]-peaktimes[0])
propvel = L/(peaktimes[Nidx-1]-peaktimes[0])
outfile = open(outfilename_vel,'w')
propvel_SI = propvel*1e-3
print('Propagation velocity (L/(peaktimes[Nidx-1]-peaktimes[0])):', propvel_SI)
outfile.write('%.16e' % propvel_SI)
outfile.close()