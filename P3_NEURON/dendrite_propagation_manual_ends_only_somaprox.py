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
testmodel = 478513437#  488462965#478513407#     480633479#  496497595 #  
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
    v_init = -86.5 # mV
elif testmodel==480633479:
    cm_soma = 0.704866 # 0.704866118957
    cm_dend = 0.704866 # 0.704866118957
    cm_axon = 0.704866 # 0.704866118957
    v_init = -96.8
elif testmodel==478513407: # Should have I=0.17
    cm_soma = 1.0
    cm_dend = 1.0
    cm_axon = 1.0
    v_init = -83.7
elif testmodel==478513437:
    cm_soma = 2.34539964752
    cm_dend = 2.34539964752
    cm_axon = 2.34539964752
    v_init = -86.8

# You can touch these instead:
cm_soma = 2.0
#cm_dend = cm_soma #5.0
#cm_axon = cm_soma #3.00603


Ls     = []
propvels = [] # List because we may encounter infs

lengthfilename = "Allen_test_changecapacitance/figures/%i/endlengths.txt" % testmodel
lengthfile = open(lengthfilename,'r') # Reading in all the info
lines = lengthfile.readlines()
for line in lines:
    words = line.split()
    Ls.append(float(words[0]))
lengthfile.close()
Nidx = len(Ls)
Nseg = Nidx#+1 # Bug
Nsegm = Nseg  # Lazy fix

folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/dendritepropagation/' # Yeah, this works for all.
outfilename_avg = folder+'idur%i_iamp' % idur+str(iamp)+'_cmspr' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_somaprox_avgpropvel_ends.txt'
outfilename_all = folder+'idur%i_iamp' % idur+str(iamp)+'_cmspr' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_somaprox_allpropvels_ends.txt'
figname_vel_dends = folder+'idur%i_iamp' % idur+str(iamp)+'_cmspr' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_somaprox_alldendpropvels_ends.png'
outfile_all = open(outfilename_all,'w')
outfilename_avg_Vmax = folder+'idur%i_iamp' % idur+str(iamp)+'_cmspr' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_somaprox_avgpropvel_Vmax_ends.txt'
outfilename_all_Vmax = folder+'idur%i_iamp' % idur+str(iamp)+'_cmspr' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_somaprox_allpropvels_Vmax_ends.txt'

filename = folder+"idur%i_iamp" % idur + str(iamp)+"_cmspr" + str(cm_soma) + "_cmd" + str(cm_dend) + "_cma"+ str(cm_axon) + "_vinit"+str(v_init)+"_wRa_V_somaandend.txt"
file  = open(filename,'r')
lines = file.readlines()
Nt    = len(lines)
times = np.zeros(Nt)
V     = np.zeros((Nseg,Nt))
peaktimes  = np.zeros(Nseg)
max_values = np.zeros(Nseg)

i=0
infnumbers  = 0
for line in lines: # Looping over dendrites
    words = line.split()
    if len(words)!=0: # words[0]:time, words[1]:V, words[2]:[Ca^2+]
        times[i]= float(words[0])
        V[0,i]  = float(words[1])
        for j in range(2,Nseg+1): # Looping over segments of dendrites
            V[j-1,i]  = float(words[j])
    i+=1
file.close()

for i in range(Nseg):        
    #peaktime = -100
    # Ugh, find the time when V is max instead...
    #for i in range (1,len(V)-1):
    #    if V[i-1]<V[i] and V[i+1]<V[i] and V:
    #        peaktimes[j] = times[i]
    #        break
    # WARNING! WARNING! NB! OOPS! OBS!:
    Vthis = V[i,:]
    max_value = max(Vthis) # Can be dangerous if I have more than one peak. But I'll make sure I don't.
    max_values[i] = max_value
    peaktimes[i]  = times[Vthis.argmax()]#times[Vthis.index(max_value)]#.index is for list only
        
# File and plot names
figname = folder+'idur%i_iamp' % idur+str(iamp)+'_cmspr' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vinit'+str(v_init)+ '_wRa_V_somaprox_enddend_prop.png'
figname_vel = folder+'idur%i_iamp' % idur+str(iamp)+'_cmspr' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_V_somaprox_enddend_propvel.png'
# Outfilename:
outfilename_vel = folder+'idur%i_iamp' % idur+str(iamp)+'_cmspr' + str(cm_soma) + '_cmd' + str(cm_dend) + '_cma'+ str(cm_axon) +'_vi'+str(v_init)+ '_wRa_V_somaprox_enddend_propvel.txt'
        
# Plotting:
    
plt.figure(figsize=(6,5))
plt.plot(peaktimes)
plt.xlabel('Index')
plt.ylabel('Peak time [ms]')
plt.title('Peak time for signal on dendrite ends')
plt.tight_layout()
plt.savefig(figname)

# Calculate speed & write all propvels together:
vel = np.zeros(Nidx)
for i in range (1,Nsegm):
    propvel_SI = Ls[i-1]/(peaktimes[i]-peaktimes[0])*1e-3 # To get it in m/s
    vel[i-1]   = propvel_SI
    if math.isinf(propvel_SI)==False:
        propvels.append(propvel_SI)
        outfile_all.write('%.16e\n' % propvel_SI)
    else:
        infnumbers+=1
outfile_all.close()

    
plt.figure(figsize=(6,5))
plt.plot(vel,'o')
plt.xlabel('Dendrite ends')
plt.ylabel('Propagation velocity [$10^{-3}$m/s]')
plt.title('Propagation velocity for signal in each point along dendrites')
plt.tight_layout()
plt.savefig(figname_vel)
    
# Double check:
print('Peaktimes:', peaktimes)

# Average and rms
propvel_avg, propvel_rms = avg_and_rms(propvels)
outfile_avg = open(outfilename_avg,'w')
# Main results:
outfile_avg.write('%.16f %.16f\n' % (propvel_avg, propvel_rms))
# Metadata:
outfile_avg.write('Number of dendrites: %i\n' % Nidx)
outfile_avg.write('propvel=inf: %i out of %i results\n' % (infnumbers,Nidx))
outfile_avg.close()

plt.figure(figsize=(6,5))
plt.plot(propvels,'o')
plt.xlabel('Dendrite index')
plt.ylabel('Propagation velocity [m/s]')
plt.title('Propagation velocity for signal in each point along different dendrites')
plt.tight_layout()
plt.savefig(figname_vel_dends)
print('cm_soma:',cm_soma)

plt.figure(figsize=(6,5))
for i in range(Nseg):
    plt.plot(times,V[i,:])
plt.xlabel('Time [ms]')
plt.ylabel('Membrane potential [mV]')
plt.title('Time vs membrane potential')
plt.tight_layout()
#plt.show()

print('V[0,:]',V[0,:])
print('V[-1,:]',V[-1,:])


### Vmax:
outfile_avg_Vmax = open(outfilename_avg_Vmax,'w')
outfile_all_Vmax = open(outfilename_all_Vmax,'w')

maxval_dend = max_values[1:]
Vmax_avg, Vmax_rms = avg_and_rms(maxval_dend)
outfile_avg_Vmax.write('%.8e %.8e' % (Vmax_avg, Vmax_rms))
outfile_avg_Vmax.close()

for i in range(len(maxval_dend)):
    outfile_all_Vmax.write('%.8e' % maxval_dend[i])
outfile_all_Vmax.close()
# Percentage of Vsoma?