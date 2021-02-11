from pylab import *
import matplotlib.pyplot as plt  # To plot
import numpy as np
import math
import random                    # Important

# Loop over hitprob? Convert to spacing? (Can always do later)
hitprob = 0.9 # Probability of hitting obstacle (and be sent back)

Nsteps    = 100   # Increase #Maybe for other values of hitprob
Nreal     = 100000
hitprev   = False # Should the system have a longer memory? # If I have only one, I don't actually need this? ###
dist_sq   = np.zeros(Nsteps+1)
positions = np.zeros(Nsteps+1) # What to do with all of these? # Average. Need to set it as an array
times     = np.arange(Nsteps+1)
steps     = np.zeros(Nsteps+1)    # probably stupid for very large Ns # Do I need this now?
hit_times = []      # These do not make sense. Maybe as arrays. Maybe I won't bother.
hit_times_pos = []  #

filename_D = 'hitprob'+str(hitprob)+'_Nsteps%i_Nreal%i_D.txt' %(Nsteps,Nreal)
filename_msd = 'hitprob'+str(hitprob)+'_Nsteps%i_Nreal%i_msd.txt' %(Nsteps,Nreal) # Should I bother to write the raw data to file? It's not like I have a lot of data point.

file_D   = open(filename_D,'w')
file_msd = open(filename_msd,'w')

for j in range(Nreal):
    step      = 0.5
    position  = 0
    for i in range(Nsteps):
        moveprob = random.random()
        if hitprev==True: # It hit something last time step, so it has gotten a little kick and will not move back to the obstacle immediately. # Possibly
            hitprev=False
            #print('hitprev true. i=',i+1, '; step before=',step)
            step*=random.random() # Does not seem to work
            #print('Step after=',step)
        else: # Move normally
            if moveprob>hitprob:
                step = 2*random.random()-1 # Return float between -1 and 1. Open for using randint instead (yields -1 and 1)
            else: # Hits something NOW
                hitprev=True
                step=-step
                #print('HIT! i=',i+1, '; step=',step)
                #hit_times.append(i)
                #hit_times_pos.append(position) # Redundant, but probably not a big deal
        #print('i=',i+1,'; step:', step)
        steps[i+1]+=step
        position+=step
        positions[i+1]+=position
        dist_sq[i+1]+=position*position

# Does this work for arrays?:
steps/=Nreal
positions/=Nreal
dist_sq/=Nreal

plt.figure(figsize=(6,5))
plt.plot(times,positions)
#plt.plot(hit_times, hit_times_pos,'.')
plt.xlabel(r't')
plt.ylabel(r'x')
plt.title(r'Position vs time')
plt.tight_layout()
plt.show()
#plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.plot(times,dist_sq)
plt.xlabel(r't')
plt.ylabel(r'MSD')
plt.title(r'MSD vs time')
plt.tight_layout()
plt.show()
#plt.savefig(plotname)

# Fit from 6 onwards actually seems ok. Still, use 10 to be safe
coeffs, covs = polyfit(times[15:], dist_sq[15:], 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a = coeffs[0]
b = coeffs[1]
D = a/2.
rms_D = np.sqrt(covs[0,0])/2.
rms_b = np.sqrt(covs[1,1])
    

file_D.write('%.16f %.16f' % (D,rms_D))
file_D.close()

for i in range(Nsteps+1):
    file_msd.write('%i %.16f\n' % (times[i],dist_sq[i]))
file_msd.close()
