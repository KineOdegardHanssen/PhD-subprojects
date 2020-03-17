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
if printevery==True:
    diffusion_end = '_D.txt'
else:
    diffusion_end = '_printevery%i_D.txt' % printevery
end_end       = totwalk_end+diffusion_end

types     = ['randomwalk', 'hardpotwalk', 'hardpotmc', 'potential', 'nearwall_rw', 'nearwall_mc', 'randomwalk1D']
filenames = [folder+'randomwalk_'+end_end, folder+'hardpotwalk_'+end_end, folder+'hardpotmc_'+end_end, folder+'pot_sigma%.3f_exp%.3f_factor%.3f_beta%.1f_' % (sigma, power, factor, beta)+end_end, 'nearwallrw_R2_Nsteps%i_Nreal%i_maxstartdist%i_printevery%i_D.txt' % (Nsteps, Nreal, maxstartdist, printevery), 'nearwallmc_R2_Nsteps%i_Nreal%i_maxstartdist%i_printevery%i_D.txt' % (Nsteps, Nreal, maxstartdist, printevery), 'rw1D_R2_Nsteps%i_Nreal%i_printevery%i_D.txt' % (Nsteps, Nreal, printevery)]
Nfiles      = len(filenames)
dateandtime = []

Ds = np.zeros(Nfiles)
kinks = []

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
    
    # Extract the time when the file was made:
    timeline   = lines[-3] # This should be the line with the time step
    #print('timeline:', timeline)
    dateandtime.append(timeline)
    
    # Extract info on kinks:
    thesekinks = []
    kinkline   = lines[-1] 
    kinks.append(kinkline)
    #kinklist   = kinkline.split()
    #for j in range(1,len(kinklist)): # Skip the first word, because that is a string
    #    thesekinks.append(int(kinklist[j]))    
    #kinks.append(kinklist)
    infile.close()
    
outfilename = folder+'/Compare/Ds_pot_sigma%.3f_exp%.3f_factor%.3f_beta%.1f_maxstartdist%i_Nsteps%i_Nreal%i' % (sigma, power, factor, beta, maxstartdist, Nsteps, Nreal)
outfile     = open(outfilename, 'w') 

for i in range(Nfiles):
    outfile.write('%s %.16f\n' % (types[i], Ds[i]))
outfile.write('\nFiles used for reading:\n')
for i in range(Nfiles):
    outfile.write('----------------------------------------------------\n')
    outfile.write('%s' % filenames[i])
    outfile.write('\n')
    outfile.write(dateandtime[i])
    outfile.write(kinks[i])
    outfile.write('\n')
    #thesekinks = kinks[i]
    #for j in range(len(thesekinks)):
    #   outfile.write(' %i' % thesekinks[j])
    #outfile.write('\n')
outfile.close()


 
    
