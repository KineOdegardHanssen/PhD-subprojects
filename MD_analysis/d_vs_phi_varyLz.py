import matplotlib.pyplot as plt                     # To plot
from pylab import *
import numpy as np

Nchains = 9
Nbeads  = 100
hsub    = 1
hcyl    = 100
r       = 0.5
psigma  = 1

# Find Lz
Lzs = []
d = []
inlocation  = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' + str(psigma) + '/Nocut/'
intfilename = inlocation+'maxzs_vs_d.txt'
infile      = open(infilename,'r')

lines = infile.readlines()
for line in lines:
    words = line.split()
    if len(words)>1: # Should always be the case, but I do not know about the last line
        d.append(float(words[0]))
        Lzs.append(float(words[1]))
d        = np.array(d)
spacings = np.array(spacings)

def give_phiblock_wsubstr(d,Lz):
    return 1 - (pi*r**2*hcyl+d**2*hsub)/(d**2*Lz)

def give_phibead_wsubstr(d,Lz):
    return 1 - (4./3*pi*r**3*(Nbeads+0.5*d**2))/(d**2*Lz)

def give_phiblock_chainsonly(d,Lz):
    return 1 - (pi*r**2*hcyl)/(d**2*Lz)

def give_phibead_chainsonly(d,Lz):
    return 1 - (4./3*pi*r**3*Nbeads)/(d**2*Lz)


### Plotting:
phi_block_ws = give_phiblock_wsubstr(d,Lzs)
phi_bead_ws  = give_phibead_wsubstr(d,Lzs)
phi_block_co = give_phiblock_chainsonly(d,Lzs)
phi_bead_co  = give_phibead_chainsonly(d,Lzs)
filename_ws = 'd_vs_phi_wsubstrate_varyLz.txt'
file_ws = open(filename_ws, 'w')
file_ws.write('d phi(blockmodel) phi(beadmodel)\n')
for i in range(len(d)):
    file_ws.write('%.2f %.16f %.16f\n' % (d[i], phi_block_ws[i], phi_bead_ws[i]))
file_ws.close()

filename_wos = 'd_vs_phi_withoutsubstrate_varyLz.txt'
file_wos = open(filename_wos, 'w')
file_wos.write('d phi(blockmodel) phi(beadmodel)\n')
for i in range(len(d)):
    file_wos.write('%.2f %.16f %.16f\n' % (d[i], phi_block_co[i], phi_bead_co[i]))
file_wos.close()