from skimage import morphology, measure, io, util   # For performing morphology operations (should maybe import more from skimage...)
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt               # To plot
from scipy.optimize import curve_fit
from pathlib import Path
import numpy as np
import random
import math
import time
import os
import glob

def halving_distances(Nx):
    Ny = Nx; Nz = Nx
    voxmat     = np.ones((Nx,Ny,Nz))
    polewidth  = int(64/4) # 64=2^6, good for halving and not too big
    voxelstart = 0
    voxelend   = polewidth-1 # To account for the change in indices
    while polewidth>1:
        voxelend = int(voxelstart+polewidth) # This should work...
        Nelem    = polewidth-1
        print('\nslab start:', voxelstart)
        print('slab end:', voxelend)
        print('slab width:', polewidth)
        voxmat[voxelstart:voxelend,:,:] = 0
        voxelstart = int(voxelend+polewidth) # Because
        polewidth = int(polewidth/2)
    if polewidth==1: # To avoid indexation problems
        voxmat[voxelstart,:,:] = 0
        voxelend = int(voxelstart+polewidth-1)
        print('\nslab start:', voxelstart)
        print('slab end:', voxelend)
        print('slab width:', polewidth)
    voxmat[Nx-1,:,:] = 0 # To 'close' the system so that we don't have to wonder/worry about the BCs
    return voxmat

def equal_distances(Nx,blockwidth): 
    Ny = Nx; Nz = Nx
    voxmat     = np.ones((Nx,Ny,Nz))
    voxelstart = 0
    voxelend   = blockwidth # To account for the change in indices
    while voxelend<Nx:
        print('\nslab start:', voxelstart)
        print('slab end:', voxelend)
        print('slab width:', blockwidth)
        voxmat[voxelstart:voxelend,:,:] = 0
        voxelstart = int(voxelend+blockwidth)
        voxelend = int(voxelstart+blockwidth) # This should work...
        print('next start:', voxelstart)
    voxmat[Nx-1,:,:] = 0 # To 'close' the system so that we don't have to wonder/worry about the BCs
    return voxmat

def spherepore(Nx,radius): 
    Ny = Nx; Nz = Nx
    voxmat      = np.zeros((Nx,Ny,Nz))
    ball        = morphology.ball(radius)
    #inverseball = morphology.ball(radius)
    #ball        = util.invert(inverseball)
    xb          = 2*radius+1 # How the ball is defined
    center      = int(Nx/2) 
    print('ball:', ball)
    voxmat[center:center+xb,center:center+xb,center:center+xb] = ball
    print('max(voxmat):', max(voxmat[center+radius, center+radius,:]))
    print('voxmat[center-radius, center-radius,:]:',voxmat[center-radius, center-radius,:])
    print('Elements in ball:', np.sum(np.sum(np.sum(ball))))
    print('Elements in voxmat:', np.sum(np.sum(np.sum(voxmat))))
    return voxmat

start_time = time.process_time()

### Grid setting
# Now, how I do choose the dimensions of the matrix...? In testing (at least for the slabs), I will see traces of the lenghts of the sheet dimensions. Maybe I need to divide by the cross-section area.
Nvoxels = 61#63
Nx      = Nvoxels
Ny      = Nvoxels
Nz      = Nvoxels
voxN    = Nx*Ny*Nz
radius  = 3
equaldistance = 6
textdist      = 'three'#'six'

isdisthalving = False
isdistequal   = False
isspherepore  = True
# Foldername
if isdisthalving:
    outfilename_base  = 'halving_distances_vox_matrix_timestep0'#'fractal_by_halving'
    foldername        = 'Voxelmatrices/Frames_halving_distances/'
    voxmat            = halving_distances(Nx)
if isdistequal:
    outfilename_base  = 'equal_distances_%s_vox_matrix_timestep0'  % textdist#'fractal_by_halving'
    foldername        = 'Voxelmatrices/Frames_equal_distances_%s/' % textdist
    voxmat            = equal_distances(Nx,equaldistance)
if isspherepore:
    outfilename_base  = 'spherepore_%s_vox_matrix_timestep0'  % textdist#'fractal_by_halving'
    foldername        = 'Voxelmatrices/Frames_spherepore_%s/' % textdist
    voxmat            = spherepore(Nx,radius)

outarray = np.zeros((voxN,4))


### Make the folder for the file if it does not already exist
p = Path(foldername)
p.mkdir(exist_ok=True)

counter = 0
prev    = 0
for i in range(Nvoxels):
    counter += 1
    value    = voxmat[i,0,0]
    if value!=prev:
        counter = 1
        print('--------------------------------------')
    print('Index', i, ', value:', value  , ' elements in group:', counter)
    prev    = value
print('Done with loop')
#print('voxmat[37:37,0,0]:',voxmat[37:37,0,0])

### Plotting?
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
pos = np.where(voxmat==1)
#print('pos set!')
ax.scatter(pos[0], pos[1], pos[2], c='black')
#print('Should plot any minute now')
#plt.savefig(plotname)
#print('Have plotted')
plt.show()

print('Have plotted')
### Opening files
outfilename_npy  = foldername + outfilename_base
outfilename_txt  = foldername + outfilename_base+'.txt'
print('Have set file names')

#outfile_txt  = open(outfilename_txt,'w')
np.save(outfilename_npy,voxmat)

print('Have saved to file')

new_time = time.process_time()
'''
for i in range(voxN): # Why this no work?
   outfile_txt.write('%i %.16f\n' % (i,outarray[i,3]))
outfile_txt.close()
'''
print("Script finished.")
