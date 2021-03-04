import numpy as np
import matplotlib.pyplot as plt

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

# Changing values of membrane capacitance:
cm = 3

# Change current:
idur = 100 # 100 # 1000 #  ms # 1
iamp = -0.5 #  0.264  # 0.41 # -0.5 # 1  # nA # 0.5 #
Ra   = 150

folder = 'Results/%i/IStim/current_idur%i_iamp'% (testmodel,idur)+str(iamp)+'/'
infilename = folder+'basaa_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit'%Ra+str(v_init)+'_V.txt'
infile = open(infilename,'r')
lines = infile.readlines()
N = len(lines) # Do I need this?

time = []
V    = []

for line in lines:
    words = line.split()
    time.append(float(words[0]))
    V.append(float(words[1]))

plt.figure(figsize=(6,5))
plt.plot(time,V)
plt.xlabel(r'$t$ [ms]')
plt.ylabel(r'$V$ [mV]')
plt.title(r'$V$ vs $t$, cm=%.1f' %cm)
plt.tight_layout()
plt.show()