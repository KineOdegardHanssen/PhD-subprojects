from pylab import *
import matplotlib.pyplot as plt  # To plot
import numpy as np
import math
import random                    # Important

r_hardsphere = 0.8    # Not entirely sure what this should be
rh2 = r_hardsphere**2 # Should be a little bit faster to work with squared values

# Setting the system. Nx=Ny=3.
N = 3; Nx = N; Ny = N; Ntarget = Nx*Ny # In case I want to extend my script later
d = 3 # Spacing. Also determines the system size
xmin = -Nx/2*d
xmax = Nx/2*d
ymin = -Ny/2*d
ymax = Ny/2*d
Lx = xmax-xmin
Ly = ymax-ymin
# Points to avoid:
targets = np.zeros((Ntarget,2)) # 9 targets with x- and y-coordinate
dstartx = -(Nx-1)/2.*d
dstarty = -(Ny-1)/2.*d
targetcounter = 0
for i in range(Nx):
    for j in range(Ny):
        targets[targetcounter,0]=dstartx+i*d
        targets[targetcounter,1]=dstarty+j*d
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
        distvec = thispos-targets[j,:]
        dist2 = np.dot(distvec,distvec)
        if dist2<rh2:
            hit_times.append(i)
            # Perform reflection # Ugh, what about the PBCs? # And ugh, what about the crossings?
            dist1 = np.sqrt(dist2)
            print('dist1:',dist1, '; r:',r)
            dist_reflected = abs(r_hardsphere-dist1) # Will move a bit before it hits the obst. Do not need to hit center-on, but this is probably good enough.
            frac_refl = 2*dist_reflected#/r
            thisx -=1.2*stepx #= positions[i-1,0] # 
            thisy -=1.2*stepy #= positions[i-1,1] #
            # Check for PBCs again:
            if thisx>xmax:
                thisx-=Lx
                print('In hit. Crossing, subtr. Lx. thisx:',thisx, '; frac_refl*stepx:',frac_refl*stepx)
            elif thisx<xmin:
                thisx+=Lx
                print('In hit. Crossing, adding Lx. thisx:',thisx, '; frac_refl*stepx:',frac_refl*stepx)
            if thisy>ymax:
                thisy-=Ly
                print('In hit. Crossing, subtr. Ly. thisy:',thisy, '; frac_refl*stepy:',frac_refl*stepy)
            elif thisy<ymin:
                thisy+=Ly
                print('In hit. Crossing, adding Ly. thisy:',thisy, '; frac_refl*stepy:',frac_refl*stepy)
            break # Cannot hit more than one obstacle in a round
    positions[i,0]=thisx
    positions[i,1]=thisy
    #print('i=',i+1,'; step:', step)
    steplens.append(r)
    dist_this = np.linalg.norm(positions[i,:]-position_0)
    dist[i]= dist_this
    dist_sq[i]= dist_this*dist_this
    times.append(i)

plt.figure(figsize=(6,5))
plt.plot(times,dist)
plt.xlabel(r't')
plt.ylabel(r'Distance from start')
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

fig, ax = plt.subplots()
#plt.figure(figsize=(6,5))
plt.plot(positions[:,0],positions[:,1], '-o')
# The alpha parameter makes the color transparent
for i in range(Ntarget):
    circle = plt.Circle((targets[i,0], targets[i,1]), r_hardsphere, color='r', alpha=0.7)
    ax.add_artist(circle) #???
plt.xlabel('x')
plt.ylabel('y')
plt.title('Trajectory')
plt.tight_layout()
plt.show()

print('xmin:',xmin)
print('xmax:',xmax)
print('ymin:',ymin)
print('ymax:',ymax)
