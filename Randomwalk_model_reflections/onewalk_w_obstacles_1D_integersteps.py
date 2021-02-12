from pylab import *
import matplotlib.pyplot as plt  # To plot
import numpy as np
import math
import random                    # Important

# Loop over hitprob? Convert to spacing? (Can always do later)
hitprob = 0.9 # Probability of hitting obstacle (and be sent back)
print('hitprob:',hitprob)

Nsteps    = 100   # Increase #Maybe for other values of hitprob
Nreal     = 100000
hitprev   = False # Should the system have a longer memory? # If I have only one, I don't actually need this? ###
dist_sq   = np.zeros(Nsteps+1)
positions = np.zeros(Nsteps+1) # What to do with all of these? # Average. Need to set it as an array
times     = np.arange(Nsteps+1)
steps     = np.zeros(Nsteps+1)    # probably stupid for very large Ns # Do I need this now?
hit_times = []      # These do not make sense. Maybe as arrays. Maybe I won't bother.
hit_times_pos = []  #


step      = 1
position  = 0
hitprob_running = hitprob
for i in range(Nsteps):
    moveprob = random.random()
    if hitprev==True: # It hit something last time step, so it has gotten a little kick and will not move back to the obstacle immediately. # Possibly
        hitprev=False
        print('hitprev true. i=',i, '; step before=',step)
        stepprev = step
        step=2*random.randint(0,1)-1
        print('i:',i,'; stepprev:',stepprev, '; step:',step)
        hitprob_running=hitprob
        if step==-stepprev: # Going to hit the target again:
            print('in test, i:',i)
            hitprob_running=1 # HAVE to hit target at next step
            #hit_times.append(i+1)
            #hit_times_pos.append(position+step)
            print('Hit by return. i+1=',i+1, '; step=',step)
    else: # Move normally
        if moveprob>hitprob_running:
            step = 2*random.randint(0,1)-1 # Return integer -1 or 1
        else: # Hits something NOW
            hitprev=True
            step=-step
            hitprob_running=hitprob
            print('HIT! i=',i, '; step=',step)
            hit_times.append(i)
            hit_times_pos.append(position) # Redundant, but probably not a big deal
    print('hitprob_running:',hitprob_running)
    print('i+1=',i+1,'; step:', step)
    steps[i+1]+=step
    position+=step
    positions[i+1]+=position
    dist_sq[i+1]+=position*position


plt.figure(figsize=(6,5))
plt.plot(times,positions)
plt.plot(hit_times, hit_times_pos,'.')
plt.xlabel(r't')
plt.ylabel(r'x')
plt.title(r'Position vs time')
plt.tight_layout()
plt.show()

plt.figure(figsize=(6,5))
plt.plot(times,dist_sq)
plt.xlabel(r't')
plt.ylabel(r'MSD')
plt.title(r'MSD vs time')
plt.tight_layout()
plt.show()