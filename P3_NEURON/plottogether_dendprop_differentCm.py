import numpy as np
import matplotlib.pyplot as plt

dend = 6
Nsec = 2 # Safeguard

# Varycm: All, soma or dendrite
#varywhichcm = 'a' # All # Only for perisomatic so far
varywhichcm = 's' # Soma
#varywhichcm = 'd' # Dendrite

## Choose model(s) to test:
#testmodel = 515175347 #Axonabundance (Hundreds of axon sections, branched dends, at least one part of the neuron is not connected to the rest) # AXON CONNECTED TO DENDRITE!!!
#testmodel = 501282059 #Originalplayer (Not everything is connected to the rest. 156 axon sec., 35 dend. sec.) # AXON CONNECTED TO DENDRITE!!!
testmodel = 496497595 #Developmentcell (One short axon, branched dendrites, everything seems to be connected, no errors)
#testmodel = 497232392
#testmodel = 496538958
#testmodel = 497230463
#testmodel = 497233075
#testmodel = 497233271 # Better Cm's, but I worry this won't run properly
#testmodel = 488462965 # PERISOMATIC model of Developmentcell
#testmodel = 478513407 # Perisomatic version of 497233271, but only yields one peak
#testmodel = 480633479 # Perisomatic version of 497230463, only one peak
#testmodel = 478513437 # Perisomatic version of 497233075
testmodelname = 'neur_%i' % testmodel
all_models    = [testmodelname]

if testmodel==480633479:
    v_init = -96.8#-83.7#-90#-86.5# # Have a test here too
elif testmodel==496497595:
    v_init = -86.5
elif testmodel==488462965:
    v_init = -86.5 # Maybe I should have changed this...
elif testmodel==497230463:
    v_init = -90
elif testmodel==497233075:
    v_init = -90
elif testmodel==478513437:
    v_init = -86.8
elif testmodel==478513407:
    v_init = -83.7
elif testmodel==497233271:
    v_init = -90

# Defaulting to original values:
# DO NOT TOUCH THESE!
# SET THEM BELOW INSTEAD!
if testmodel==496497595:
    cm_soma = 1.14805
    cm_dend = 9.98231
    cm_axon = 3.00603
    if dend==6:
        Nsec = 17 # Starts at 0
elif testmodel==497233271:
    cm_soma = 0.783229
    cm_dend = 1.94512
    cm_axon = 8.25387
elif testmodel==497230463:
    cm_soma = 1.23729
    cm_dend = 2.57923
    cm_axon = 5.02697
elif testmodel==497233075:
    cm_soma = 1.64168
    cm_dend = 2.83035
    cm_axon = 9.98442
elif testmodel==488462965:
    cm_soma = 3.31732779736 # Strange values, right?
    cm_dend = 3.31732779736
    cm_axon = 3.31732779736
elif testmodel==478513407:
    cm_soma = 1.0
    cm_dend = 1.0
    cm_axon = 1.0
elif testmodel==480633479:
    cm_soma = 0.704866118957
    cm_dend = 0.704866118957
    cm_axon = 0.704866118957
elif testmodel==478513437:
    cm_soma = 2.34539964752
    cm_dend = 2.34539964752
    cm_axon = 2.34539964752


# Changing values of membrane capacitance:
cms = [0.01,3.0]#[1.14805,3.0] # Include option 'default?'

# Change current:
idur = 2.0#100 # 100 # 1000 #  ms # 1
iamp = 1.0#-0.5 #  0.264  # 0.41 # -0.5 # 1  # nA # 0.5 #

### Reading in and plotting #######
folder = 'figures/%i/current_idur%i_iamp'% (testmodel,idur)+str(iamp)+'/dendritepropagation/'

fig=plt.figure(figsize=(19,9))

### First capacitance: #######
if varywhichcm=='a':
   cm_soma = cms[0]
   cm_dend = cms[0]
   cm_axon = cms[0]
   cmsize  = str(cms[0])
   textsnippet = 'whole neuron'
   plotname = folder+'idur%i_iamp'%idur+str(iamp)+'_varyallcms'+str(cms[0])+'and'+str(cms[1])+'_vinit'+str(v_init)+'_wRa_d%i.png' % dend
elif varywhichcm=='s':
   cm_soma = cms[0]
   cmsize  = str(cms[0])
   textsnippet = 'soma'
   plotname = folder+'idur%i_iamp'%idur+str(iamp)+'_varycmsoma'+str(cms[0])+'and'+str(cms[1])+'_cmd'+str(cm_dend)+'_cma'+str(cm_axon)+'_vinit'+str(v_init)+'_wRa_V_d%i.png' % dend
elif varywhichcm=='d':
   cm_dend = cms[0]
   cmsize  = str(cms[0])
   textsnippet = 'dendrite'
   plotname = folder+'idur%i_iamp'%idur+str(iamp)+'_varycmdend'+str(cms[0])+'and'+str(cms[1])+'_cms'+str(cm_soma)+'_cma'+str(cm_axon)+'_vinit'+str(v_init)+'_wRa_V_d%i.png' % dend

fig.add_subplot(121, xlabel="time (ms)", ylabel="Membrane potential (mV)", title=r"Potential vs time for %s capcacitance $C_m=$%s" % (textsnippet,cmsize))

# Sub-plot
for i in range(Nsec):
    infilename = folder+'idur%i_iamp'%idur+str(iamp)+'_cms'+str(cm_soma)+'_cmd'+str(cm_dend)+'_cma'+str(cm_axon)+'_vinit'+str(v_init)+'_wRa_V_dsec%i_d%i.txt' % (i,dend)
    infile = open(infilename,'r')
    lines = infile.readlines()
    N = len(lines) # Do I need this?
    
    time = []
    V    = []
    
    for line in lines:
        words = line.split()
        time.append(float(words[0]))
        V.append(float(words[1]))
    plt.plot(time,V,label='Section %i' %i)
plt.legend(loc='upper right')

### Second capacitance: #######
if varywhichcm=='a':
   cm_soma = cms[1]
   cm_dend = cms[1]
   cm_axon = cms[1]
   cmsize  = str(cms[1])
elif varywhichcm=='s':
   cm_soma = cms[1]
   cmsize  = str(cms[1])
   textsnippet = 'soma'
elif varywhichcm=='d':
   cm_dend = cms[1]
   cmsize  = str(cms[1])
   textsnippet = 'dendrite'
fig.add_subplot(122, xlabel="time (ms)", ylabel="Membrane potential (mV)", title=r"Potential vs time for %s capcacitance $C_m=$%s" % (textsnippet,cmsize))

# Sub-plot
for i in range(Nsec):
    infilename = folder+'idur%i_iamp'%idur+str(iamp)+'_cms'+str(cm_soma)+'_cmd'+str(cm_dend)+'_cma'+str(cm_axon)+'_vinit'+str(v_init)+'_wRa_V_dsec%i_d%i.txt' % (i,dend)
    infile = open(infilename,'r')
    lines = infile.readlines()
    N = len(lines) # Do I need this?
    
    time = []
    V    = []
    
    for line in lines:
        words = line.split()
        time.append(float(words[0]))
        V.append(float(words[1]))
    plt.plot(time,V,label='Section %i' %i)
plt.legend(loc='upper right')

fig.savefig(plotname)