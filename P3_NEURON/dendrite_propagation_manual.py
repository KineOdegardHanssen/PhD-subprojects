import matplotlib.pyplot as plt
import numpy

##### Adjustable parameters/specifications #########
# File/simulation selection:
testmodel = 496497595
idur = 1 # ms
iamp = 1.0 # nA
v_init = -86.5 # mV

dendnr = 6
idxs = [37,38,39,40,41,42,43,44,45,46]
Nidx = len(idxs)
####################################################

# Defaulting to original values:
# DO NOT TOUCH THESE!
cm_soma = 1.14805
cm_dend = 9.98231
cm_axon = 3.00603

# You can touch these instead:
#cm_soma = 1.14805
#cm_dend = 9.98231
#cm_axon = 3.00603

peaktimes = numpy.zeros(Nidx)


folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/dendritepropagation/'

for j in range(Nidx):
    filename = folder+'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) +'_vinit'+str(v_init)+ '_addedRa_V_dendsec%i_dend%i.txt' % (j,dendnr)
    
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
    for i in range (1,len(V)-1):
        if V[i-1]<V[i] and V[i+1]<V[i]:
            peaktimes[j] = times[i]
            break


# Plot name:
figname = folder+'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) +'_vinit'+str(v_init)+ '_addedRa_V_dend%i_propagation.png' % dendnr

# Plotting:

plt.figure(figsize=(6,5))
plt.plot(idxs, peaktimes)
plt.xlabel('Index')
plt.ylabel('Peak time [ms]')
plt.title('Peak time for signal along dendrite %i' % dendnr)
plt.tight_layout()
plt.savefig(figname)
