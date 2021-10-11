import numpy as np
import matplotlib.pyplot as plt
import math

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
v_init = -70 

# Changing values of membrane capacitance:
cm = 5.0

# Change current:
idur = 1000#100 # 100 # 1000 #  ms # 1
iamp = -10.0 #  0.264  # 0.41 # -0.5 # 1  # nA # 0.5 #
Ra   = 100

folder = 'Results/Istim/current_idur%i_iamp'% idur+str(iamp)+'/'
infilename = folder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit'%Ra+str(v_init)+'_V.txt'
infile = open(infilename,'r')
lines = infile.readlines()
N = len(lines) # Do I need this?

time = []
V    = []

for line in lines:
    words = line.split()
    time.append(float(words[0]))
    V.append(float(words[1]))

Vmax = max(V)
Vmin = min(V)
dV = Vmax-Vmin
dV1e = dV*math.exp(-1)
V1e  = Vmin+dV1e
print('Vmax:',Vmax)
print('Vmin:',Vmin)
print('V1e:',V1e)


for i in range(len(V)):
    if V[i]-V1e<0:
        print('V:',V[i-1],'i-1:',i-1,'; t:', time[i-1])
        print('V:',V[i],'i:',i,'; t:', time[i])
        plot_t = time[i]
        break


plt.figure(figsize=(6,5))
plt.plot(time,V)
plt.plot(plot_t,V1e,'o')
plt.xlabel(r'$t$ [ms]')
plt.xlabel(r'$V$ [mV]')
plt.title(r'$V$ vs $t$')
plt.tight_layout()
plt.show()