from pylab import *
import matplotlib.pyplot as plt  # To plot
import numpy as np
import math
import random                    # Important

rh = 0.8    # Size of obstacle. Not entirely sure what this should be. rh=0.65?

# Setting the system. Nx=Ny=3.
N = 3; Nx = N; Ny = N; Ntarget = Nx # In case I want to extend my script later
d = 3 # Spacing. Also determines the system size
xmin = -Nx/2*d
xmax = Nx/2*d
ymin = -Ny/2*d
ymax = Ny/2*d
Lx = xmax-xmin
Ly = ymax-ymin
# Points to avoid:
targets = np.zeros(Ntarget) # 3 targets with x-coordinate
dstartx = -(Nx-1)/2.*d
dstarty = -(Ny-1)/2.*d
targetcounter = 0
for i in range(Nx):
    targets[targetcounter]=dstartx+i*d
    targetcounter+=1

# Should I operate with velocities? Anders said so.
# I should average over many walks. Will study ONE first.

Nsteps    = 100   # Increase
Nreal     = 1     # Increase (not yet implemented)
hitprev   = False # Should the system have a longer memory? # If I have only one, I don't actually need this? ###
r         = 0.5     # Default # Remove? ###
step      = [0.5,0]
thisx     = 1
thisy     = 1
positions = np.zeros((Nsteps+1,2))
position_0 = np.array([1,1])
positions[0,0] = 1
positions[0,1] = 1
dist           = np.zeros(Nsteps+1)
dist_sq        = np.zeros(Nsteps+1)
times          = [0]
steplens       = []    # probably stupid for very large Ns
angles         = []    # probably stupid for very large Ns
hit_times      = []

z0 = position_0[0]

for i in range(1,Nsteps+1):
    # Draw random numbers
    r = random.random() # Return float between 0 and 1.
    a = random.uniform(0,2*np.pi)
    stepx = r*np.cos(a)
    stepy = r*np.sin(a)
    # Temporary update of positions:
    thisx = positions[i-1,0]+stepx
    thisy = positions[i-1,1]+stepy
    # Then apply the PBCs: # Should this be in a function?
    if thisx>xmax:
        thisx-=Lx
        print('Crossing, subtr. Lx. thisx:',thisx)
    elif thisx<xmin:
        thisx+=Lx
        print('Crossing, adding Lx. thisx:',thisx)
    if thisy>ymax:
        thisy-=Ly
        print('Crossing, subtr. Ly. thisy:',thisy)
    elif thisy<ymin:
        thisy+=Ly
        print('Crossing, adding Ly. thisy:',thisy)
    thispos = np.array([thisx,thisy])
    # Check if it hits a target
    for j in range(Ntarget):
        xpos = thispos[0]-targets[j]
        sep  = abs(xpos)
        if sep<rh:
            # PROBABILITY OF PERFORMING REFLECTION?
            hit_times.append(i)
            # Perform reflection # Ugh, what about the PBCs? # And ugh, what about the crossings?
            thisx = positions[i-1,0]
            thisy = positions[i-1,1]+stepy # Assuming that it hits the obstable midways in the time step.
            # Check for PBCs again:
            if thisx>xmax:
                thisx-=Lx
                print('In hit. Crossing, subtr. Lx. thisx:',thisx)
            elif thisx<xmin:
                thisx+=Lx
                print('In hit. Crossing, adding Lx. thisx:',thisx)
            if thisy>ymax:
                thisy-=Ly
                print('In hit. Crossing, subtr. Ly. thisy:',thisy)
            elif thisy<ymin:
                thisy+=Ly
                print('In hit. Crossing, adding Ly. thisy:',thisy)
            break # Cannot hit more than one obstacle in a round
    positions[i,0]=thisx
    positions[i,1]=thisy
    #print('i=',i+1,'; step:', step)
    steplens.append(r)
    dist_this = positions[i,1]-z0
    dist[i]= dist_this
    dist_sq[i]= dist_this*dist_this
    times.append(i)

plt.figure(figsize=(6,5))
plt.plot(times,dist)
plt.xlabel(r't')
plt.ylabel(r'Distance from start, z-direction')
plt.title(r'Distance vs time')
plt.tight_layout()
plt.show()
#plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.plot(times,positions[:,0], '-o')
plt.xlabel(r't')
plt.ylabel(r'x')
plt.title(r'x vs time')
plt.tight_layout()
plt.show()

plt.figure(figsize=(6,5))
plt.plot(times,positions[:,1], '-o')
plt.xlabel(r't')
plt.ylabel(r'z')
plt.title(r'z vs time')
plt.tight_layout()
plt.show()

plt.figure(figsize=(6,5))
plt.plot(times,dist_sq)
plt.xlabel(r't')
plt.ylabel(r'Distance$^2$ from start (in z-direction)')
plt.title(r'Distance$^2$ vs time')
plt.tight_layout()
plt.show()


fig, ax = plt.subplots()
#plt.figure(figsize=(6,5))
plt.plot(positions[:,0],positions[:,1], '-o')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Trajectory')
plt.tight_layout()
plt.show()

for target in targets:
    print('target:',target)

print('xmin:',xmin)
print('xmax:',xmax)
print('zmin:',ymin)
print('zmax:',ymax)
