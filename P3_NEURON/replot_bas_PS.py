import numpy as np
import matplotlib.pyplot as plt

## Choose model(s) to test:
#####testmodel = 480633479 # Perisomatic version of 497230463. Weird model.
#testmodel = 488462965 # PERISOMATIC model of Developmentcell
#testmodel = 478513407 # Perisomatic version of 497233271, but only yields one peak
testmodel = 478513437 # Perisomatic version of 497233075
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
    v_init = -86.5 # Wrong name... -86.8
elif testmodel==478513407:
    v_init = -83.7
elif testmodel==497233271:
    v_init = -90

# Changing values of membrane capacitance:
cm = 0.5

# Change current:
idur = 1000 # 100 # 1000 #  ms # 1
iamp = 0.5 #  # nA # 0.5 #
Ra   = 150

#folder = 'C:/Users/Kine/Documents/Projects_PhD/P3_PNN_Capacitance/Cellemodeller/Ball-and-stick models/BAS_somaPV_dendPV/Results/%i/IStim/current_idur%i_iamp'% (testmodel,idur)+str(iamp)+'/'
folder = 'Results/%i/IStim/current_idur%i_iamp'% (testmodel,idur)+str(iamp)+'/'
infilename = folder+'basps_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit'%Ra+str(v_init)+'_V.txt'
infile = open(infilename,'r')
lines = infile.readlines()
N = len(lines) # Do I need this?

time = []
V    = []

for line in lines:
    words = line.split()
    time.append(float(words[0]))
    V.append(float(words[1]))

peak_heights = []
vmax = max(V) 
vmin = min(V) 
deltav = vmax-vmin
vthr  = vmax-0.15*deltav # If there is a peak above this value, we count it
vprev = vthr-40 # A peak never kicks in at initiation, even if I change vthr
Npeaks = 0
for i in range (1,len(V)-1):  
    if V[i-1]<V[i] and V[i+1]<V[i] and V[i]>vthr:
        Npeaks+=1
        peak_heights.append(V[i])
avg_peak_height = np.mean(peak_heights)
rms = 0
for i in range(Npeaks):
    rms += (peak_heights[i]-avg_peak_height)**2
rms = np.sqrt(rms/(Npeaks-1))
print(Npeaks, ' peaks for model ', testmodel, ', current ', iamp)
print('vmax:', vmax)
print('peak heights:',peak_heights)
print('avg peak height:',avg_peak_height,'; rms:', rms)
print('testmodel:',testmodel)

plt.figure(figsize=(6,5))
plt.plot(time,V)
plt.xlabel(r'$t$ [ms]')
plt.ylabel(r'$V$ [mV]')
plt.title(r'$V$ vs $t$, cm=%.1f' %cm)
plt.tight_layout()
plt.show()