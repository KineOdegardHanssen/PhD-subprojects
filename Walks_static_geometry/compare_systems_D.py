from mpl_toolkits.mplot3d import Axes3D             # Plotting in 3D
import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob
import copy

# Choosing file
Nsteps = 20000
Nreal  = 1000
Nsect  = 5
beta   = 3
intd   = 1
intsigma = 10
Nblocks  = 3
sigma    = intsigma
power    = 6
factor   = 1
printevery   = 100
maxstartdist = 10

printall    = False ### PRINTALL REQUIRES SOME CHANGE!!!

# Save fig or show fig
savefig = False

### Setting file names
folder        = 'PBC/sigma%i_d%i/Nblocks%i/' % (intsigma, intd, Nblocks)
totwalk_end   = 'R2_Nsteps%i_Nreal%i' % (Nsteps, Nreal)
diffusion_end = '_D.txt'
end_end       = totwalk_end+diffusion_end

types     = ['randomwalk', 'hardpotwalk', 'hardpotmc', 'potential', 'nearwall_rw', 'nearwall_mc', 'randomwalk1D']
filenames = [folder+'randomwalk_'+end_end, folder+'hardpotwalk_'+end_end, folder+'hardpotmc_'+end_end, folder+'pot_sigma%.3f_exp%.3f_factor%.3f_beta%.1f_' % (sigma, power, factor, beta)+end_end, 'nearwall_rw_R2_Nsteps%i_Nreal%i_maxstartdist%i_printevery%i_D.txt' % (Nsteps, Nreal, maxstartdist, printevery), 'nearwall_mc_R2_Nsteps%i_Nreal%i_maxstartdist%i_printevery%i_D.txt' % (Nsteps, Nreal, maxstartdist, printevery), 'rw1D_R2_Nsteps%i_Nreal%i_printevery%i_D.txt' % (Nsteps, Nreal, printevery)]
Nfiles = len(filenames)

Ds = np.zeros(Nfiles)

for i in range(Nfiles):
    ### Read file totwalk
    filename = filenames[i]
    infile   = open(filename, 'r')
    
    lines = infile.readlines()
    Nlines = len(lines)
    
    # Extract second to last line and use that as D:
    firstline = lines[0]
    NDs       = int(firstline.split()[1])
    line = lines[NDs] # Extracting D from the last fit. The first line does not contain any Ds, so the Ds start at line number 1 (and not 0). Hence, this is the correct line.
    D    = float(line.split()[2])
    Ds[i] = D
    infile.close()
    
outfilename = folder+'/Compare/Ds_pot_sigma%.3f_exp%.3f_factor%.3f_beta%.1f_maxstartdist%i_Nsteps%i_Nreal%i' % (sigma, power, factor, beta, maxstartdist, Nsteps, Nreal)
outfile     = open(outfilename, 'w') 

for i in range(NFiles):
    outfile.write('%s %.16f\n' % (types[i], Ds[i]))
outfile.close()


 
    
