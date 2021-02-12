from pylab import *
import matplotlib.pyplot as plt  # To plot
import numpy as np
import math
import random                    # Important

# Plot as a function of hitprob?
hitprob = 0.2 # Probability of hitting obstacle (and be sent back)

# Threshold for colliding with the same obstacle again:
collthr = -0.9

# Should I operate with velocities? Anders said so.
# I should average over many walks. Will study ONE first.

Nsteps    = 100   # Increase
Nreal     = 1     # Increase (not yet implemented)
hitprev   = False # Should the system have a longer memory? # If I have only one, I don't actually need this? ###
step      = [1,0]
position  = [0,0]
positions = np.zeros((Nsteps+1,2))
dist      = np.zeros(Nsteps+1)
dist_sq   = np.zeros(Nsteps+1)
times     = [0]
angles    = []    # probably stupid for very large Ns
hit_times = []
hit_times_dist = []
hit_times_x    = []
hit_times_y    = []

hitprob_running = hitprob
for i in range(Nsteps):
    moveprob = random.random()
    if hitprev==True: # It hit something last time step, so it has gotten a little kick and will not move back to the obstacle immediately. # Possibly
        hitprev=False
        hitprob_running = hitprob
        stepprev = step
        a = random.uniform(0,2*np.pi)
        step = np.array([np.cos(a),np.sin(a)])
        print('hitprev true. i=',i+1, '; step before=',stepprev, '; step now:',step)
        dotprod = np.dot(step,stepprev)
        if dotprod<collthr:
            hitprob_running = 1
            print('Rehit, i+1=',i+1,'; dotprod:',dotprod)
        print('Step after=',step)
    else: # Move normally
        if moveprob>hitprob_running:
            a = random.uniform(0,2*np.pi)
            step = np.array([np.cos(a),np.sin(a)])
        else: # Hits something NOW
            hitprev=True
            # This here?: hitprob_running = hitprob
            # How to determine how it hits # Randomize that too?
            hitdir = random.randint(0,2) # OR should I draw the angle of the surface orientation?
            if hitdir==0:            # Is reflection just as simple?
                step=np.array([-step[0],-step[1]])
            elif hitdir==1:
                step = np.array([-step[0],step[1]])
            elif hitdir==2:
                step = np.array([step[0],-step[1]])
            print('HIT! i=',i+1, '; step=',step)
            hit_times.append(i)
            hit_times_dist.append(dist[i]) # Redundant, but probably not a big deal
            hit_times_x.append(position[0])
            hit_times_y.append(position[1])
    print('i=',i+1,'; step:', step)
    position+=step
    positions[i+1,:]=position
    dist_this = np.linalg.norm(position)
    dist[i+1]= dist_this
    dist_sq[i+1]= dist_this*dist_this
    times.append(i+1)

plt.figure(figsize=(6,5))
plt.plot(times,dist)
plt.plot(hit_times, hit_times_dist,'.')
plt.xlabel(r't')
plt.ylabel(r'Distance from start')
plt.title(r'Distance vs time')
plt.tight_layout()
plt.show()
#plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.plot(times,positions[:,0])
plt.plot(hit_times, hit_times_x,'.')
plt.xlabel(r't')
plt.ylabel(r'x')
plt.title(r'x vs time')
plt.tight_layout()
plt.show()

plt.figure(figsize=(6,5))
plt.plot(times,positions[:,1])
plt.plot(hit_times, hit_times_y,'.')
plt.xlabel(r't')
plt.ylabel(r'y')
plt.title(r'y vs time')
plt.tight_layout()
plt.show()


plt.figure(figsize=(6,5))
plt.plot(times,dist_sq)
plt.xlabel(r't')
plt.ylabel(r'Distance$^2$ from start')
plt.title(r'Distance$^2$ vs time')
plt.tight_layout()
plt.show()
