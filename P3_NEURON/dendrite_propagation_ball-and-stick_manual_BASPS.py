import matplotlib.pyplot as plt
import numpy as np

# Plotting vs cm too
testmodel = 478513407#478513437#488462965 #
##### Adjustable parameters/specifications #########
# File/simulation selection:
idur = 2  # ms
iamp = 1.0 # nA
v_init = -65 # mV
L      = 980 # 1000 with the soma, but will exclude that in the calculations
smallCms = True

Ra   = 150
if testmodel==488462965:
    Ra = 29
elif testmodel==478513407:
    Ra = 100
elif testmodel==478513437:
    Ra = 352

if smallCms==False:
    cms = [0.5,1,2,3,4,4.5,5,5.5,6,7,8,9,10]
else: 
    cms = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.6,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
idxs    = np.array([1,13,25,37,49]) # Dendrite only
Ncm     = len(cms)
Nidx    = len(idxs)
dL      = L/float(Nidx)
peaktimes = np.zeros(Nidx)
propvels_SI = np.zeros(Ncm)
####################################################

folder = 'Results/%i/IStim/current_idur'%testmodel+str(idur)+'_iamp'+str(iamp)+'/dendritepropagation/'
filename_cms = folder+'basps_varycm_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_propvel_vs_cm'
figname_vel_vs_cm = folder+'basps_varycm_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_propvel_vs_cm'

if smallCms==True:
    filename_cms      = filename_cms+'_smallCm.txt'
    figname_vel_vs_cm = figname_vel_vs_cm+'_smallCm.png'
else:
    filename_cms      = filename_cms+'.txt'
    figname_vel_vs_cm = figname_vel_vs_cm+'.png'

file_cms = open(filename_cms,'w')
for i in range(Ncm):
    cm = cms[i]
    for j in range(1,Nidx):
        filename = folder+'basps_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_V_seg%i.txt' %j # For easy access
        
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
        
        max_value = max(V)
        peaktimes[j] = times[V.index(max_value)]
    
    # File and plot names
    
    figname = folder+'basps_cm'+str(cm)+'_idur'+str(idur)+'_iamp'+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_V_dend_propagation.png'
    figname_vel = folder+'basps_cm'+str(cm)+'_idur'+str(idur)+'_iamp'+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_V_dend_propvel.png'
    # Outfilename:
    outfilename_vel = folder+'basps_cm'+str(cm)+'_idur'+str(idur)+'_iamp'+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_V_dend_propvel.txt'
    
     
    # Plotting:
    
    plt.figure(figsize=(6,5))
    plt.plot(idxs, peaktimes)
    plt.xlabel('Index')
    plt.ylabel('Peak time [ms]')
    plt.title('Peak time for signal along dendrite')
    plt.tight_layout()
    plt.savefig(figname)
    
    # Calculate speed:
    vel = np.zeros(Nidx-1)
    for k in range (Nidx-2):
        vel[k] = dL/(peaktimes[k+1]-peaktimes[k])
    
    plt.figure(figsize=(6,5))
    plt.plot(idxs[:Nidx-1], vel)
    plt.xlabel('Index')
    plt.ylabel('Propagation velocity [$10^{-3}$m/s]')
    plt.title('Propagation velocity for signal in each point along dendrite')
    plt.tight_layout()
    plt.savefig(figname_vel)
    
    # Double check:
    print('Peaktimes:', peaktimes)
    
    # Other approach: print L/(peaktimes[Nidx-1]-peaktimes[0])
    propvel = L/(peaktimes[Nidx-1]-peaktimes[0])
    outfile = open(outfilename_vel,'w')
    propvel_SI = propvel*1e-3
    propvels_SI[i] = propvel_SI
    print('Propagation velocity (L/(peaktimes[Nidx-1]-peaktimes[0])):', propvel_SI)
    outfile.write('%.16e' % propvel_SI)
    outfile.close()
    file_cms.write('%.2f %.16e\n' % (cm,propvel_SI))
file_cms.close()

# Plot propvel vs cm 
plt.figure(figsize=(6,5))
plt.plot(cms, propvels_SI,'-o')
plt.xlabel(r'Cell capacitance $C_m$')
plt.ylabel('Signal velocity [$10^{-3}$m/s]')
plt.title('Signal velocity along dendrite vs $C_m$ ')
plt.tight_layout()
plt.savefig(figname_vel_vs_cm)