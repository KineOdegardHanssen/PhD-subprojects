from pylab import *
import matplotlib.pyplot as plt  # To plot
import numpy as np
import math
import random                    # Important

r_hardsphere = 1.0   # Not entirely sure what this should be
rh2 = r_hardsphere**2 # Should be a little bit faster to work with squared values
reflfac = 1.2
maxlen  = 0.1

# Setting the system. Nx=Ny=3.
N = 3; Nx = N; Ny = N; Nz = N; Ntarget = Nx*Ny*Nz # In case I want to extend my script later
d = 25 # Spacing. Also determines the system size
xmin = -Nx/2*d
xmax = Nx/2*d
ymin = -Ny/2*d
ymax = Ny/2*d
zmin = -Nz/2*r_hardsphere
zmax = Nz/2*r_hardsphere
Lx = xmax-xmin
Ly = ymax-ymin
Lz = zmax-zmin
# Points to avoid:
targets = np.zeros((Ntarget,3)) # 9 targets with x- and y-coordinate
dstartx = -(Nx-1)/2.*d
dstarty = -(Ny-1)/2.*d
dstartz = -(Nz-1)/2.*r_hardsphere
targetcounter = 0
for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz):
            targets[targetcounter,0]=dstartx+i*d
            targets[targetcounter,1]=dstarty+j*d
            targets[targetcounter,2]=dstartz+k*r_hardsphere
            targetcounter+=1


print('d:',d)
# Should I operate with velocities? Anders said so.
# I should average over many walks. Will study ONE first.

Nsteps    = 1000 # Increase?
Nreal     = 2#10000 # Increase?
hitprev   = False # Should the system have a memory?
# Total displacement
dist           = np.zeros(Nsteps+1)
dist_sq        = np.zeros(Nsteps+1)
dist_uw        = np.zeros(Nsteps+1)
dist_sq_uw     = np.zeros(Nsteps+1)
#Parallel displacement
dist_par       = np.zeros(Nsteps+1)
dist_par_sq    = np.zeros(Nsteps+1)
dist_par_uw    = np.zeros(Nsteps+1)
dist_par_sq_uw = np.zeros(Nsteps+1)
#Perpendicular displacement
dist_ort       = np.zeros(Nsteps+1)
dist_ort_sq    = np.zeros(Nsteps+1)
dist_ort_uw    = np.zeros(Nsteps+1)
dist_ort_sq_uw = np.zeros(Nsteps+1)
# Time array
times          = np.linspace(0,Nsteps,Nsteps+1)
steplens       = []    # probably stupid for very large Ns
angles         = []    # probably stupid for very large Ns
hit_times      = []

for k in range(Nreal):
    print('k:',k)
    thisx     = 1
    thisy     = 1
    positions = np.zeros((Nsteps+1,3))
    position_0 = np.array([1,1,1])
    positions[0,0] = 1
    positions[0,1] = 1
    positions[0,2] = 1
    positions_uw = np.zeros((Nsteps+1,3))
    positions_uw[0,0] = 1
    positions_uw[0,1] = 1
    positions_uw[0,2] = 1
    for i in range(1,Nsteps+1):
        # Draw random numbers
        u = random.random()*maxlen**3 # Return float between 0 and maxlen**3.
        r = u**(1./3)                 # Return float between 0 and maxlen.
        v = 2*random.random()-1       # Return float between -1 and 1
        vroot = np.sqrt(1-v*v)
        theta = random.random()*2*np.pi # Limits: [). Supposed to be ()? Float between 0 and 2pi.
        stepx = r*vroot*np.cos(theta)
        stepy = r*vroot*np.sin(theta)
        stepz = r*v
        # Temporary update of positions:
        thisx = positions[i-1,0]+stepx
        thisy = positions[i-1,1]+stepy
        thisz = positions[i-1,2]+stepz
        # Then apply the PBCs: # Should this be in a function?
        if thisx>xmax:
            thisx-=Lx
        elif thisx<xmin:
            thisx+=Lx
        if thisy>ymax:
            thisy-=Ly
        elif thisz<zmin:
            thisz+=Lz
        elif thisz<zmin:
            thisz+=Lz
        thispos = np.array([thisx,thisy,thisz])
        # Check if it hits a target
        for j in range(Ntarget):
            distvec = thispos-targets[j,:]
            dist2 = np.dot(distvec,distvec)
            if dist2<rh2:
                hit_times.append(i)
                # Do not move:
                thisx=positions[i-1,0]
                thisy=positions[i-1,1]
                thisz=positions[i-1,2]
                stepx=0
                stepy=0
                stepz=0
                break # Cannot hit more than one obstacle in a round
        positions[i,0]=thisx
        positions[i,1]=thisy
        positions[i,2]=thisz
        positions_uw[i,0]=positions_uw[i-1,0]+stepx
        positions_uw[i,1]=positions_uw[i-1,1]+stepy
        positions_uw[i,2]=positions_uw[i-1,2]+stepz
        #print('i=',i+1,'; step:', step)
        dist_this        = np.linalg.norm(positions[i,:]-position_0)
        dist_uw_this     = np.linalg.norm(positions_uw[i,:]-position_0)
        dist_par_this    = np.linalg.norm(positions[i,0:1]-position_0[0:1])
        dist_par_uw_this = np.linalg.norm(positions_uw[i,0:1]-position_0[0:1])
        dist_ort_this    = abs(positions[i,2]-position_0[2])
        dist_ort_uw_this = abs(positions_uw[i,2]-position_0[2])
        ### Into array:
        #      Total
        dist[i]+= dist_this
        dist_sq[i]+= dist_this*dist_this
        dist_uw[i]+= dist_uw_this
        dist_sq_uw[i]+= dist_uw_this*dist_uw_this
        #      Parallel
        dist_par[i]+= dist_par_this
        dist_par_sq[i]+= dist_par_this*dist_par_this
        dist_par_uw[i]+= dist_par_uw_this
        dist_par_sq_uw[i]+= dist_par_uw_this*dist_par_uw_this
        #      Orthogonal
        dist_ort[i]+= dist_ort_this
        dist_ort_sq[i]+= dist_ort_this*dist_ort_this
        dist_ort_uw[i]+= dist_ort_uw_this
        dist_ort_sq_uw[i]+= dist_ort_uw_this*dist_ort_uw_this
#Total
dist/=Nreal
dist_sq/=Nreal
dist_uw/=Nreal
dist_sq_uw/=Nreal
#Parallel
dist_par/=Nreal
dist_par_sq/=Nreal
dist_par_uw/=Nreal
dist_par_sq_uw/=Nreal
#Orthogonal
dist_ort/=Nreal
dist_ort_sq/=Nreal
dist_ort_uw/=Nreal
dist_ort_sq_uw/=Nreal

plt.figure(figsize=(6,5))
plt.plot(times,dist_uw)
plt.xlabel(r't')
plt.ylabel(r'Distance from start')
plt.title(r'Distance vs time')
plt.tight_layout()
plt.show()


## Fits
# Total
coeffs, covs = polyfit(times[10:], dist_sq_uw[10:], 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a = coeffs[0]
b = coeffs[1]
D = a/6.
rms_D = np.sqrt(covs[0,0])/6.
rms_b = np.sqrt(covs[1,1])
# Parallel
coeffs, covs = polyfit(times[10:], dist_par_sq_uw[10:], 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a = coeffs[0]
b = coeffs[1]
Dpar = a/4.
rms_Dpar = np.sqrt(covs[0,0])/4.
rms_bpar = np.sqrt(covs[1,1])
# Orthogonal
coeffs, covs = polyfit(times[10:], dist_ort_sq_uw[10:], 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a = coeffs[0]
b = coeffs[1]
Dort = a/2.
rms_Dort = np.sqrt(covs[0,0])/2.
rms_bort = np.sqrt(covs[1,1])


filename_D = '3D_d'+str(d)+'_rsphere'+str(r_hardsphere)+'_Nsteps%i_Nreal%i_norefl_maxlen'%(Nsteps,Nreal)+str(maxlen)+'_D.txt'


filename_msd = '3D_d'+str(d)+'_rsphere'+str(r_hardsphere)+'_Nsteps%i_Nreal%i_norefl_maxlen'%(Nsteps,Nreal)+str()+'_msd.txt'

plotname_msd = '3D_d'+str(d)+'_rsphere'+str(r_hardsphere)+'_Nsteps%i_Nreal%i_norefl_maxlen'%(Nsteps,Nreal)+str()+'_msd.png'

file_D   = open(filename_D,'w')
file_msd = open(filename_msd,'w')

file_D.write('%.16f %.16f %.16f %.16f %.16f %.16f' % (D,rms_D, Dpar, rms_Dpar, Dort, rms_Dort))
file_D.close()

for i in range(Nsteps+1):
    file_msd.write('%i %.16f %.16f %.16f\n' % (times[i],dist_sq_uw[i],dist_par_sq_uw[i],dist_ort_sq_uw[i],))
file_msd.close()

print('positions_uw[:,0]:',positions_uw[:,0])
print('positions_uw[:,1]:',positions_uw[:,1])

print('min(positions_uw[:,0]):',min(positions_uw[:,0]))
print('min(positions_uw[:,1]):',min(positions_uw[:,1]))
print('max(positions_uw[:,0]):',max(positions_uw[:,0]))
print('max(positions_uw[:,1]):',max(positions_uw[:,1]))


plt.figure(figsize=(6,5))
plt.plot(times,dist_sq_uw, label='Total') # Should I include rms-values? Get big...
plt.plot(times,dist_par_sq_uw, label='Parallel')
plt.plot(times,dist_ort_sq_uw, label='Orthogonal')
plt.xlabel(r't')
plt.ylabel(r'Distance$^2$ from start')
plt.title(r'Distance$^2$ vs time')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(plotname_msd)
plt.show()

print('targets:',targets)
print('len(targets):',len(targets))
print('zmin:',zmin)
print('dstartz:',dstartz)

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
plt.title('Trajectory, unwrapped')
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
plt.title('Trajectory, wrapped')
plt.tight_layout()
plt.show()
'''

