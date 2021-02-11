from pylab import *
import matplotlib.pyplot as plt  # To plot
import numpy as np
import math
import random                    # Important

# Plot as a function of hitprob?
hitprob = 0.2 # Probability of hitting obstacle (and be sent back)

# Should I operate with velocities? Anders said so.
# I should average over many walks. Will study ONE first.

Nsteps    = 100   # Increase
Nreal     = 1     # Increase (not yet implemented)
hitprev   = False # Should the system have a longer memory? # If I have only one, I don't actually need this? ###
dirprev   = 1     # Default # Remove? ###
step      = 0.5
position  = 0
positions = [0]
times     = [0]
steps     = []    # probably stupid for very large Ns
hit_times = []
hit_times_pos = []

for i in range(Nsteps):
    moveprob = random.random()
    if hitprev==True: # It hit something last time step, so it has gotten a little kick and will not move back to the obstacle immediately. # Possibly
        hitprev=False
        print('hitprev true. i=',i+1, '; step before=',step)
        step*=random.random() # Does not seem to work
        print('Step after=',step)
    else: # Move normally
        if moveprob>hitprob:
            step = 2*random.random()-1 # Return float between -1 and 1. Open for using randint instead (yields -1 and 1)
        else: # Hits something NOW
            hitprev=True
            step=-step
            print('HIT! i=',i+1, '; step=',step)
            hit_times.append(i)
            hit_times_pos.append(position) # Redundant, but probably not a big deal
    print('i=',i+1,'; step:', step)
    steps.append(step)
    position+=step
    positions.append(position)
    times.append(i+1)

plt.figure(figsize=(6,5))
plt.plot(times,positions)
plt.plot(hit_times, hit_times_pos,'.')
plt.xlabel(r't')
plt.ylabel(r'x')
plt.title(r'Position vs time')
plt.tight_layout()
plt.show()
#plt.savefig(plotname)
    

    

