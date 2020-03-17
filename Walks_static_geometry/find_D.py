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

randomwalk  = False
hardpot_rw  = False
hardpot_mc  = False
potential   = False
nearwall_rw = False
nearwall_mc = False
onedim      = True
printall    = False

# Save fig or show fig
savefig = True

# Kinks?
kinks_in = [1200, 2300, 3300] # nearwall_rw #[2000,3500] # potential, Nsteps20000, Nreal1000, Nsect5, sigma10, d=1
kinks    = []
for i in range(len(kinks_in)):
    thiskink = int(floor(kinks_in[i]/printevery))
    kinks.append(thiskink)

### Setting file names
prefix = ''
if randomwalk==True:
    prefix = 'randomwalk_'
elif hardpot_rw==True:
    prefix = 'hardpotwalk_'
elif hardpot_mc==True:
    prefix = 'hardpotmc_'
elif potential==True:
    prefix = 'pot_sigma%.3f_exp%.3f_factor%.3f_beta%.1f_' % (sigma, power, factor, beta) # I set it to three decimals. Maybe a bit much, but perhaps it is better to be careful.

folder = 'PBC/sigma%i_d%i/Nblocks%i/' % (intsigma, intd, Nblocks)
infilename_totwalk  = folder + prefix + 'R2_Nsteps%i_Nreal%i' % (Nsteps, Nreal)
infilename_sections = folder + prefix + 'R2_Nsteps%i_Nreal%i_Npart%i' % (Nsteps, Nreal, Nsect)
outfilename          = infilename_totwalk + '_D.txt'

#infilename_totwalk  = 'Nsteps20000_Nreal1000_sigma2_d10_Nblocks3_randomwalk_PBC_R2_basic_mc' #'Nsteps20000_Nreal100_Nblocks3_intsigma2_intd10_pot_exp6.000000_sigma1.000000_factor1.000000_PBC_R2_basic_mc' #'Nsteps20000_Nreal100_hardpot_R2_basic_mc'
#infilename_sections = 'Nsteps20000_Nreal1000_sigma2_d10_Nblocks3_Npart5_randomwalk_PBC_R2_basic_mc' #'Nsteps20000_Nreal100_Nblocks3_intsigma2_intd10_Npart5_pot_exp6.000000_sigma1.000000_factor1.000000_PBC_R2_basic_mc' #'Nsteps20000_Nreal100_Npart5_hardpot_R2_basic_mc'

plotname_sections_wfit = infilename_sections + '_step_vs_R2_wfit'

if printall==False:
    plotname_sections_wfit = plotname_sections_wfit + '_printevery%i.png' % printevery
    infilename_sections    = infilename_sections    + '_printevery%i' % printevery
    infilename_totwalk     = infilename_totwalk     + '_printevery%i' % printevery
    outfilename            = infilename_totwalk     + '_D.txt'
else:
    plotname_sections_wfit = plotname_sections_wfit + '.png'

if nearwall_rw==True or nearwall_mc==True or onedim==True: # Totally different file convention for these ones.
    difftype = ''
    if nearwall_rw==True:
        difftype = 'rw'
    else:
        difftype = 'mc'
    infilename_totwalk     = 'nearwall'+difftype + '_R2_Nsteps%i_Nreal%i_maxstartdist%i_printevery%i' % (Nsteps, Nreal, maxstartdist, printevery)
    infilename_sections    = 'nearwall'+difftype + '_R2_Nsteps%i_Nreal%i_maxstartdist%i_Npart%i_printevery%i' % (Nsteps, Nreal, maxstartdist, Nsect, printevery)
    if onedim==True:
        infilename_totwalk  = 'rw1D_R2_Nsteps%i_Nreal%i_printevery%i' % (Nsteps, Nreal, printevery)
        infilename_sections = 'rw1D_R2_Nsteps%i_Nreal%i_Npart%i_printevery%i' % (Nsteps, Nreal, Nsect, printevery)
    plotname_sections_wfit = infilename_sections + '.png'
    outfilename            = infilename_totwalk + '_D.txt'


### Read file totwalk
infile_totwalk = open(infilename_totwalk, 'r')

lines = infile_totwalk.readlines()
#Nlines = len(lines)
#print("Nlines:", Nlines)

# Ready lists
step    = [0] # Nothing has happened in the first time step
R2      = [0]
R2_stdv = [0]

for line in lines:
    words = line.split()
    if len(words)>2: # This should always be true, but why not make sure
        step.append(int(words[0]))
        R2.append(float(words[1]))
        R2_stdv.append(float(words[2]))
infile_totwalk.close()

step    = np.array(step)
R2      = np.array(R2)
R2_stdv = np.array(R2_stdv)


### Read file sections
infile_sections = open(infilename_sections, 'r')

lines = infile_sections.readlines()
#Nlines = len(lines)
#print("Nlines:", Nlines)

# Ready lists
step_all    = [] 
R2_all      = []
R2_stdv_all = []
Nsect       = 0  # We should know this from the file name, but we can also do it like this

for line in lines:
    words = line.split()
    if len(words)==4:
        Nsect += 1
        if int(words[2])==0: # If we start # Could use pop, but...
            step_section    = []
            R2_section      = []
            R2_stdv_section = []
        else:
            step_all.append(np.array(step_section))
            R2_all.append(np.array(R2_section))
            R2_stdv_all.append(np.array(R2_stdv_section))
            step_section    = []
            R2_section      = []
            R2_stdv_section = []
    if len(words)==3: # This should always be true, but why not make sure
        step_section.append(int(words[0]))
        R2_section.append(float(words[1]))
        R2_stdv_section.append(float(words[2]))
step_all.append(step_section)
R2_all.append(R2_section)
R2_stdv_all.append(R2_stdv_section)

infile_sections.close()

## Done with read-in

# Append all data
step_all_flat = []
R2_all_flat   = []
for i in range(Nsect):
    step_this = step_all[i]
    R2_this   = R2_all[i]
    for j in range(len(step_this)):
        step_all_flat.append(step_this[j])
        R2_all_flat.append(R2_this[j])
# Polyfit
Ds = []
a_array = []
b_array = []

# Fit of all the data
#popt, pcov = np.polyfit(step_all[0],R2_all, 1) # Not so sure if this is the right approach after all... Gives one a and b for each data set, it seems
popt = np.polyfit(step_all_flat,R2_all_flat, 1)#, cov='unscaled')
a = popt[0]
b = popt[1]
D = a/4.
fit = a*step_all[0]+b
Ds.append(D)
a_array.append(a)
b_array.append(b)

# Fit of parts
start = 0
thesteps = step_all[0] # The steps are the same for all
for i in range(len(kinks)+1):
    if i!=len(kinks):
        end = kinks[i]
    else:
        #print('In else')
        start = kinks[i-1]
        end   = int(floor(thesteps[-1]/printevery))
    length = end-start
    steps_thisfit = []
    R2s_thisfit   = []
    for j in range(Nsect):
        R2_these = R2_all[j]
        for k in range(end-start):
            #print('i =', i,' last step:', len(thesteps-1),'; start:', start, ', end:', end, 'start+k:', start+k)
            steps_thisfit.append(thesteps[start+k])
            R2s_thisfit.append(R2_these[start+k])
    popt = np.polyfit(steps_thisfit,R2s_thisfit, 1)#, cov='unscaled')
    a = popt[0]
    b = popt[1]
    D = a/4.
    #fit = a*step_all[0]+b # Should save the parameters?
    Ds.append(D)
    a_array.append(a)
    b_array.append(b)
    start = end


plt.figure(figsize=(6,5))
for i in range(Nsect):
    plt.errorbar(step_all[i], R2_all[i], yerr=R2_stdv_all[i], fmt="none", capsize=2, label='Section %i' % i)
plt.plot(step_all[0], fit, label='Line fit, all')
for i in range(len(kinks)):
    a = a_array[i+1]
    b = b_array[i+1]
    thisfit = a*step_all[0]+b
    plt.plot(step_all[0], thisfit, label='Line fit %i' % (i+1))
plt.xlabel(r'Step', fontsize=16)
plt.ylabel(r'<$R^2$>', fontsize=16)
plt.tight_layout(pad=3.0)#, w_pad=0.0, h_pad=0.5)
plt.legend(loc="lower right")
plt.title(r'<$R^2$> vs step', fontsize=16) # (last chain)', fontsize=16)
if savefig==True:
    plt.savefig(plotname_sections_wfit)
else:
    plt.show()

NDs = len(kinks)+1
outfile = open(outfilename, 'w')
outfile.write('NDs: %i' % NDs)
outfile.write('D, all: %.16f' % Ds[0])
for i in range(NDs):
    outfile.write('\nD, fit%i: %.16f' % ((i+1),Ds[i+1]))
outfile.write('\n\nKinks')
for i in range(len(kinks_in)):
    outfile.write(' %i' % kinks_in[i])
outfile.close()


