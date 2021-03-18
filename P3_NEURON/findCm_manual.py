import numpy as np
import matplotlib.pyplot as plt
import math

## Choose model(s) to test:
#testmodel = 515175347 #Axonabundance (Hundreds of axon sections, branched dends, at least one part of the neuron is not connected to the rest) # AXON CONNECTED TO DENDRITE!!!
#testmodel = 501282059 #Originalplayer (Not everything is connected to the rest. 156 axon sec., 35 dend. sec.) # AXON CONNECTED TO DENDRITE!!!
#testmodel = 496497595 #Developmentcell (One short axon, branched dendrites, everything seems to be connected, no errors)
#testmodel = 497232392
#testmodel = 496538958
#testmodel = 497230463
#testmodel = 497233075
#testmodel = 497233271 # Better Cm's, but I worry this won't run properly
testmodel = 488462965 # PERISOMATIC model of Developmentcell
#testmodel = 478513407 # Perisomatic version of 497233271, but only yields one peak
#testmodel = 480633479 # Perisomatic version of 497230463, only one peak
#testmodel = 478513437 # Perisomatic version of 497233075
testmodelname = 'neur_%i' % testmodel
all_models    = [testmodelname]

if testmodel==480633479:
    v_init = -96.8#-83.7#-90#-86.5# # Have a test here too
elif testmodel==496497595:
    A      = 2286e-12
    v_init = -86.5
elif testmodel==488462965:
    A      = 2286e-12
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
cm_soma = 0.01
cm_dend = cm_soma #2.0#
cm_axon = cm_soma #0.01

# Change current:
idur = 1000 #50 #2.0#100 # 100 # 1000 #  ms # 1
iamp = -0.5 # 1.0 # 0.264  # 0.41 # -0.5 # 1  # nA # 0.5 #
idelay = 100
testit = False # If I am unsure about the results, I can check the fit.

folder = 'figures/%i/current_idur%i_iamp'% (testmodel,idur)+str(iamp)+'/'
infilename = folder+'idur%i_iamp'%idur+str(iamp)+'_cmsoma'+str(cm_soma)+'_cmdend'+str(cm_dend)+'_cmaxon'+str(cm_axon)+'_vinit'+str(v_init)+'_addedRa.txt'
infile = open(infilename,'r')
lines = infile.readlines()
N = len(lines) # Do I need this?

time = []
V    = []
Ca   = []

for line in lines:
    words = line.split()
    time.append(float(words[0]))
    V.append(float(words[1]))
    Ca.append(float(words[2]))

Vmax = max(V)
Vmin = min(V)
dV = Vmax-Vmin
dV1e = dV*math.exp(-1)
V1e  = Vmin+dV1e
print('Vmax:',Vmax)
print('Vmin:',Vmin)
print('V1e:',V1e)


# Will perform linear interpolation near the crossing with V1e.
tbef = 0
taft = 0
Vbef = 0
Vaft = 0
for i in range(len(V)):
    if V[i]-V1e<0:
        tbef = time[i-1]
        taft = time[i]
        Vbef = V[i-1]
        Vaft = V[i]
        print('V:',V[i-1],'i-1:',i-1,'; t:', time[i-1])
        print('V:',V[i],'i:',i,'; t:', time[i])
        break

a = (Vaft-Vbef)/(taft-tbef)
t_interpolated = (V1e-Vbef+a*tbef)/a-idelay

print('tbef',tbef)
print('taft',taft)
print('t_interpolated',t_interpolated)

tau1 = tbef-idelay
tau2 = taft-idelay
tau3 = t_interpolated

print('dV+Vmax:',dV+Vmax)
if testit==True:
    N = 100*idur
    timecut = np.linspace(0,idur,N)
    fit1 = np.zeros(N)
    fit2 = np.zeros(N)
    fit3 = np.zeros(N)
    for i in range(N):
        fit1[i] = dV*math.exp(-timecut[i]/tau1)+Vmin
        fit2[i] = dV*math.exp(-timecut[i]/tau2)+Vmin
        fit3[i] = dV*math.exp(-timecut[i]/tau3)+Vmin
    timecut = np.linspace(idelay,idur+idelay,N)
    
    plt.figure(figsize=(6,5))
    plt.plot(time,V, label='Data')
    plt.plot(timecut,fit1, 'r--', label='tbef')
    plt.plot(timecut,fit2, 'g--', label='taft')
    plt.plot(timecut,fit3, 'b--', label='tinterpolated')
    plt.plot(tau1+idelay,V1e, 'ro', label='tbef')
    plt.plot(tau2+idelay,V1e, 'go', label='taft')
    plt.plot(tau3+idelay,V1e, 'bo', label='tinterpolated')
    plt.xlabel(r'$t$ [ms]')
    plt.ylabel(r'$V$ [mV]')
    plt.title(r'$V$ vs $t$')
    plt.tight_layout()
    plt.legend(loc='upper right')
    plt.show()

folder = 'figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/'
filename = folder+'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) + '_vinit'+str(v_init)+'_addedRa.txt'
outfilename = folder +'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) + '_vinit'+str(v_init)+'_addedRa_Cmtotal_manual.txt'

tau = tau3

# Conversions (to SI):
iamp*=1e-9
tau*=1e-3

dV*= 1e-3
R  = abs(dV/iamp)
C  = tau/R
Cm = C/A*100 # Conversion to get muF/cm^2

outfile = open(outfilename,'w')
outfile.write('%.12f' % Cm)
outfile.close()