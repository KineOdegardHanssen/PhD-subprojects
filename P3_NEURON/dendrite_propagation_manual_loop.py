import matplotlib.pyplot as plt
import numpy as np
import math
import sys

def avg_and_rms(x):
    N = len(x)
    avgx = np.mean(x)
    print('avgx:',avgx)
    rmsx = 0
    for i in range(N):
        rmsx += (avgx-x[i])**2
    rmsx = np.sqrt(rmsx/(N-1))
    return avgx,rmsx

##### Adjustable parameters/specifications #########
# File/simulation selection:
testmodel = 496497595 # 488462965 # 
iduraa = 2
if testmodel==496497595:
    idur = iduraa # ms
elif testmodel==488462965:
    idur = 2 # ms
iamp = 1.0 # nA
v_init = -86.5 # mV


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
#cm_soma = 0.01
cm_dend = 10.0 #cm_soma #5.0
#cm_axon = cm_soma #3.00603


Nidx   = [] # I think I will abandon using the actual segment names. Too much of a hassle and not that useful
Ls     = []
propvels = [] # List because we may encounter infs

morphfolder = 'Allen_test_changecapacitance/figures/%i/' % testmodel # this is all I need.
morphfilename = morphfolder+'%i_dendrites_section_segments.txt' % testmodel
morphfile = open(morphfilename,'r') # Reading in all the info
lines = morphfile.readlines()
for line in lines:
    words = line.split()
    Nsegm = int(words[1])
    if Nsegm==0:
        print('WARNING! EMPTY DENDRITE SECTION!')
        print('The soma or axon probably have more sections than you thought. Look into that before attempting to run this script again.')
        print('(because you will mix axon and dendrite sections now)')
        sys.exit(1) # Just exiting this mess
    Nidx.append(Nsegm)
    Ls.append(float(words[2]))
morphfile.close()

folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/dendritepropagation/' # Yeah, this works for all.
outfilename_avg = folder+'idur%i_iamp' % idur+str(iamp)+'_cms' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_avgpropvel.txt'
outfilename_all = folder+'idur%i_iamp' % idur+str(iamp)+'_cms' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_allpropvels.txt'
figname_vel_dends = folder+'idur%i_iamp' % idur+str(iamp)+'_cms' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_alldendpropvels.png'

outfile_all = open(outfilename_all,'w')
infnumbers  = 0
for k in range(len(Nidx)): # Looping over dendrites
    Nidxk = Nidx[k]
    L = Ls[k]
    dL = L/Nidxk
    peaktimes = np.zeros(Nidxk)
    for j in range(Nidxk): # Looping over segments of dendrites
        filename = folder+"idur%i_iamp" % idur + str(iamp)+"_cms" + str(cm_soma) + "_cmd" + str(cm_dend) + "_cma"+ str(cm_axon) + "_vinit"+str(v_init)+"_wRa_V_dsec%i_d%i.txt" % (j,k)
        
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
        # WARNING! WARNING! NB! OOPS! OBS!:
        max_value = max(V) # Can be dangerous if I have more than one peak. But I'll make sure I don't.
        peaktimes[j] = times[V.index(max_value)]
        
    # File and plot names
    figname = folder+'idur%i_iamp' % idur+str(iamp)+'_cms' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vinit'+str(v_init)+ '_wRa_V_dend%i_prop.png' % k
    figname_vel = folder+'idur%i_iamp' % idur+str(iamp)+'_cms' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_V_dend%i_propvel.png' % k
    # Outfilename:
    outfilename_vel = folder+'idur%i_iamp' % idur+str(iamp)+'_cms' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_V_dend%i_propvel.txt' % k
        
    # Plotting:
    
    plt.figure(figsize=(6,5))
    plt.plot(peaktimes)
    plt.xlabel('Index')
    plt.ylabel('Peak time [ms]')
    plt.title('Peak time for signal along dendrite %i' % k)
    plt.tight_layout()
    plt.savefig(figname)
    
    # Calculate speed:
    vel = np.zeros(Nidx[k]-1)
    for i in range (0,Nidx[k]-1):
        vel[i] = dL/(peaktimes[i+1]-peaktimes[i])
    
    plt.figure(figsize=(6,5))
    plt.plot(vel)
    plt.xlabel('Index')
    plt.ylabel('Propagation velocity [$10^{-3}$m/s]')
    plt.title('Propagation velocity for signal in each point along dendrite %i' % k)
    plt.tight_layout()
    plt.savefig(figname_vel)
    
    # Double check:
    print('Peaktimes:', peaktimes)
    
    # Other approach: print L/(peaktimes[Nidx-1]-peaktimes[0])
    propvel = L/(peaktimes[Nidx[k]-1]-peaktimes[0])
    outfile = open(outfilename_vel,'w')
    propvel_SI = propvel*1e-3
    print('Propagation velocity (L/(peaktimes[Nidx-1]-peaktimes[0])):', propvel_SI)
    outfile.write('%.16e' % propvel_SI)
    outfile.close()
    # Writing all propvels together:
    if math.isinf(propvel)==False:
        propvels.append(propvel_SI)
        outfile_all.write('%.16e\n' % propvel_SI)
    else:
        infnumbers+=1
outfile_all.close()

propvel_avg, propvel_rms = avg_and_rms(propvels)
outfile_avg = open(outfilename_avg,'w')
# Main results:
outfile_avg.write('%.16f %.16f\n' % (propvel_avg, propvel_rms))
# Metadata:
outfile_avg.write('Number of dendrites: %i\n' % len(Nidx))
outfile_avg.write('propvel=inf: %i out of %i results\n' % (infnumbers,len(Nidx)))
outfile_avg.close()

plt.figure(figsize=(6,5))
plt.plot(propvels,'o')
plt.xlabel('Dendrite index')
plt.ylabel('Propagation velocity [m/s]')
plt.title('Propagation velocity for signal in each point along different dendrites')
plt.tight_layout()
plt.savefig(figname_vel_dends)