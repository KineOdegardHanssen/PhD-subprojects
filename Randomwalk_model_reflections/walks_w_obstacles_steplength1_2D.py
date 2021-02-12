from pylab import *
import matplotlib.pyplot as plt  # To plot
import numpy as np
import math
import random                    # Important

# Plot as a function of hitprob?
hitprob = 0 # Probability of hitting obstacle (and be sent back)
print('hitprob:',hitprob)

# Threshold for colliding with the same obstacle again:
collthr = -0.9 # Is this too negative? # I should probably print this to the file name

Nsteps    = 100   # Increase
Nreal     = 10000  # Increase
hitprev   = False # Should the system have a longer memory? # If I have only one, I don't actually need this? ###
positions = np.zeros((Nsteps+1,2))
dist      = np.zeros(Nsteps+1)
dist_sq   = np.zeros(Nsteps+1)
times     = np.arange(0,Nsteps+1)
steplens  = []    # probably stupid for very large Ns
angles    = []    # probably stupid for very large Ns


filename_D = '2D_hitprob'+str(hitprob)+'_Nsteps%i_Nreal%i_collthr' %(Nsteps,Nreal)+str(collthr)+'_D.txt'
filename_msd = '2D_hitprob'+str(hitprob)+'_Nsteps%i_Nreal%i_collthr' %(Nsteps,Nreal)+str(collthr)+'_msd.txt' # Should I bother to write the raw data to file? It's not like I have a lot of data point.

file_D   = open(filename_D,'w')
file_msd = open(filename_msd,'w')

for j in range(Nreal):
    step     = [1,0]
    position = [0,0]
    hitprob_running = hitprob
    for i in range(Nsteps):
        moveprob = random.random()
        if hitprev==True: # It hit something last time step, so it has gotten a little kick and will not move back to the obstacle immediately. # Possibly
            hitprev=False
            hitprob_running = hitprob
            stepprev = step
            a = random.uniform(0,2*np.pi)
            step = np.array([np.cos(a),np.sin(a)])
            dotprod = np.dot(step,stepprev)
            if dotprod<collthr:
                hitprob_running = 1
        else: # Move normally
            if moveprob>hitprob_running:
                a = random.uniform(0,2*np.pi)
                step = np.array([np.cos(a),np.sin(a)])
            else: # Hits something NOW
                hitprev=True
                # How to determine how it hits # Randomize that too?
                hitdir = random.randint(0,2) # OR should I draw the angle of the surface orientation?
                if hitdir==0:            # Is reflection just as simple?
                    step=np.array([-step[0],-step[1]])
                elif hitdir==1:
                    step = np.array([-step[0],step[1]])
                elif hitdir==2:
                    step = np.array([step[0],-step[1]])
        position[0]+=step[0]
        position[1]+=step[1]
        positions[i+1,:]+=position
        dist_this = np.linalg.norm(position)
        dist[i+1]+= dist_this
        dist_sq[i+1]+= dist_this*dist_this

dist/=Nreal
dist_sq/=Nreal

plt.figure(figsize=(6,5))
plt.plot(times,dist)
plt.xlabel(r't')
plt.ylabel(r'Distance from start')
plt.title(r'Distance vs time')
plt.tight_layout()
plt.show()
#plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.plot(times,dist_sq)
plt.xlabel(r't')
plt.ylabel(r'Distance$^2$ from start')
plt.title(r'Distance$^2$ vs time')
plt.tight_layout()
plt.show()

# Fit actually seems ok for all times for d<=0.5. Need to change depending for the larger hitprob
# Strangely enough, the msd is more noisy for large times.
coeffs, covs = polyfit(times[2:], dist_sq[2:], 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a = coeffs[0]
b = coeffs[1]
D = a/4.
rms_D = np.sqrt(covs[0,0])/4.
rms_b = np.sqrt(covs[1,1])
    

file_D.write('%.16f %.16f' % (D,rms_D))
file_D.close()

for i in range(Nsteps+1):
    file_msd.write('%i %.16f\n' % (times[i],dist_sq[i]))
file_msd.close()
