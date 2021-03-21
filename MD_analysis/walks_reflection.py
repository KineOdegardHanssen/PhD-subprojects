from pylab import *
import matplotlib.pyplot as plt  # To plot
import numpy as np
import math
import random                    # Important

r_hardsphere = 0.65 #1.0    # Not entirely sure what this should be
rh2 = r_hardsphere**2 # Should be a little bit faster to work with squared values
reflfac = 1.2

# Setting the system. Nx=Ny=3.
N = 3; Nx = N; Ny = N; Ntarget = Nx*Ny # In case I want to extend my script later
d = 100 # Spacing. Also determines the system size
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

Nsteps    = 1000 # Increase?
Nreal     = 10000 # Increase?
hitprev   = False # Should the system have a memory?
dist           = np.zeros(Nsteps+1)
dist_sq        = np.zeros(Nsteps+1)
dist_uw        = np.zeros(Nsteps+1)
dist_sq_uw     = np.zeros(Nsteps+1)
times          = np.linspace(0,Nsteps,Nsteps+1)
steplens       = []    # probably stupid for very large Ns
angles         = []    # probably stupid for very large Ns
hit_times      = []

for k in range(Nreal):
    thisx     = 1
    thisy     = 1
    positions = np.zeros((Nsteps+1,2))
    position_0 = np.array([1,1])
    positions[0,0] = 1
    positions[0,1] = 1
    positions_uw = np.zeros((Nsteps+1,2))
    positions_uw[0,0] = 1
    positions_uw[0,1] = 1
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
        elif thisx<xmin:
            thisx+=Lx
        if thisy>ymax:
            thisy-=Ly
        elif thisy<ymin:
            thisy+=Ly
        thispos = np.array([thisx,thisy])
        # Check if it hits a target
        for j in range(Ntarget):
            distvec = thispos-targets[j,:]
            dist2 = np.dot(distvec,distvec)
            if dist2<rh2:
                hit_times.append(i)
                # Perform reflection # Ugh, what about the PBCs? # And ugh, what about the crossings?
                dist1 = np.sqrt(dist2)
                dist_reflected = abs(r_hardsphere-dist1) # Will move a bit before it hits the obst. Do not need to hit center-on, but this is probably good enough.
                frac_refl = 2*dist_reflected#/r
                thisx -=reflfac*stepx #= positions[i-1,0] # 
                thisy -=reflfac*stepy #= positions[i-1,1] #
                stepx =(1-reflfac)*stepx
                stepy =(1-reflfac)*stepy
                # Check for PBCs again:
                if thisx>xmax:
                    thisx-=Lx
                elif thisx<xmin:
                    thisx+=Lx
                if thisy>ymax:
                    thisy-=Ly
                elif thisy<ymin:
                    thisy+=Ly
                break # Cannot hit more than one obstacle in a round
        positions[i,0]=thisx
        positions[i,1]=thisy
        positions_uw[i,0]=positions_uw[i-1,0]+stepx
        positions_uw[i,1]=positions_uw[i-1,1]+stepy
        #print('i=',i+1,'; step:', step)
        dist_this = np.linalg.norm(positions[i,:]-position_0)
        dist_uw_this = np.linalg.norm(positions_uw[i,:]-position_0)
        #if dist_uw_this>12.7: # Should only kick in for uw coord.
        #    print('dist_uw_this:',dist_uw_this)
        dist[i]+= dist_this
        dist_sq[i]+= dist_this*dist_this
        dist_uw[i]+= dist_uw_this
        dist_sq_uw[i]+= dist_uw_this*dist_uw_this

dist/=Nreal
dist_sq/=Nreal
dist_uw/=Nreal
dist_sq_uw/=Nreal

plt.figure(figsize=(6,5))
plt.plot(times,dist_uw)
plt.xlabel(r't')
plt.ylabel(r'Distance from start')
plt.title(r'Distance vs time')
plt.tight_layout()
plt.show()
#plt.savefig(plotname)


plt.figure(figsize=(6,5))
plt.plot(times,dist_sq_uw)
plt.xlabel(r't')
plt.ylabel(r'Distance$^2$ from start')
plt.title(r'Distance$^2$ vs time')
plt.tight_layout()
plt.show()

# Fit actually seems ok for all times for d<=0.5. Need to change depending for the larger hitprob
# Strangely enough, the msd is more noisy for large times.
coeffs, covs = polyfit(times[2:], dist_sq_uw[2:], 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a = coeffs[0]
b = coeffs[1]
D = a/4.
rms_D = np.sqrt(covs[0,0])/4.
rms_b = np.sqrt(covs[1,1])


filename_D = '2D_d'+str(d)+'_rsphere'+str(r_hardsphere)+'_reflfac'+str(reflfac)+'_Nsteps%i_Nreal%i_D.txt' %(Nsteps,Nreal)
filename_msd = '2D_d'+str(d)+'_rsphere'+str(r_hardsphere)+'_reflfac'+str(reflfac)+'_Nsteps%i_Nreal%i_msd.txt' %(Nsteps,Nreal) # Should I bother to write the raw data to file? It's not like I have a lot of data point.

file_D   = open(filename_D,'w')
file_msd = open(filename_msd,'w')

file_D.write('%.16f %.16f' % (D,rms_D))
file_D.close()

for i in range(Nsteps+1):
    file_msd.write('%i %.16f\n' % (times[i],dist_sq[i]))
file_msd.close()

print('positions_uw[:,0]:',positions_uw[:,0])
print('positions_uw[:,1]:',positions_uw[:,1])

print('min(positions_uw[:,0]):',min(positions_uw[:,0]))
print('min(positions_uw[:,1]):',min(positions_uw[:,1]))
print('max(positions_uw[:,0]):',max(positions_uw[:,0]))
print('max(positions_uw[:,1]):',max(positions_uw[:,1]))

'''
fig, ax = plt.subplots()
#plt.figure(figsize=(6,5))
plt.plot(positions_uw[:,0],positions_uw[:,1], '-o')
# The alpha parameter makes the color transparent
for i in range(Ntarget):
    circle = plt.Circle((targets[i,0], targets[i,1]), r_hardsphere, color='r', alpha=0.7)
    ax.add_artist(circle) #???
plt.xlabel('x')
plt.ylabel('y')
plt.title('Trajectory')
plt.tight_layout()
plt.show()
'''