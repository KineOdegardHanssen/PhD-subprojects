import matplotlib.pyplot as plt
import numpy as np
#from scipy.optimize import curve_fit

# Function for curve_fit
#def negative_exponential(s,P): # Have prefactor A too? Maybe do this manually?
#    return np.exp(-s/P)

tol = 1e-1 # Tolerance - when V-V_last is smaller than this, cut the data

testmodel = 496497595 #
idur   = 50 # ms
idelay = 10
iamp   = -0.5 # nA
v_init = -86.5 # mV


if testmodel==496497595: # Will only implement this, probably
    cm_soma = 1.14805
    cm_dend = 9.98231
    cm_axon = 3.00603

# Changing values of membrane capacitance:
#cm_soma = 1.0#0.01
#cm_dend = cm_soma#10.0#
#cm_axon = cm_soma #0.01#

# Bruke efel til å lese inn?
# Kan også sjekke med tids-arrayen
# dt = t-t_onset
t = []
V = []

folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp'% (testmodel,idur)+str(iamp)+'/'
filename = filename = folder +'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) + '_vinit'+str(v_init)+'_addedRa.txt'
file = open(filename,'r')
lines = file.readlines()

for line in lines:
    words = line.split()
    if len(words)>0:
        t.append(float(words[0]))
        V.append(float(words[1]))
file.close()

V_last = V[-1] # Base value for testing
V_first = V[0] # Have 'equilibrated'. Do I really need this, though?

dt = [] # Delta t, not time step
V_fall = []
for i in range(len(t)):
    thist = t[i]
    thisV = V[i]
    if abs(thisV-V_last)<tol:
        break
    if thist>=idelay:
        dt.append(thist-idelay)
        V_fall.append(thisV)

dt = np.array(dt)
V_fall = np.array(V_fall)

V_fall_positive = V_fall - min(V_fall) # Subtracting the most negative number: The smallest value is 0

# For ylog:
Vfp_log = np.log(V_fall_positive)

# Should plot to verify

findex = 100#5#100  # First index
lindex = 300#40#300 # Last index
exponent = (Vfp_log[lindex]-Vfp_log[findex])/(dt[lindex]-dt[findex])
tau = -1./exponent

fit = dt*exponent+max(Vfp_log)

print(tau)
print('len(Vfp_log):',len(Vfp_log))

plt.figure(figsize=(6,5))
plt.plot(dt,Vfp_log)
plt.plot(dt,fit,'--')
plt.xlabel('Time [ms]')
plt.ylabel('log(V)') # Units make no sense here
plt.title('log(V) vs time')
plt.show()

plt.figure(figsize=(6,5))
plt.plot(t,V)
plt.xlabel('Time [ms]')
plt.ylabel('V')
plt.title('V vs time')
plt.show()

print('tau:',tau)

# Conversions (to SI):
iamp*=1e-9
tau*=1e-3
A   = 2286e-12

dV = (max(V)-min(V))*1e-3
R  = abs(dV/iamp)
C  = tau/R
Cm = C/A*100

print('R:', R)
print('C:',C)
print('Cm:',Cm)