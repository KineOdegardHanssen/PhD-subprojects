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
printevery = 100
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
savefig = False

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

#infilename_totwalk  = 'Nsteps20000_Nreal1000_sigma2_d10_Nblocks3_randomwalk_PBC_R2_basic_mc' #'Nsteps20000_Nreal100_Nblocks3_intsigma2_intd10_pot_exp6.000000_sigma1.000000_factor1.000000_PBC_R2_basic_mc' #'Nsteps20000_Nreal100_hardpot_R2_basic_mc'
#infilename_sections = 'Nsteps20000_Nreal1000_sigma2_d10_Nblocks3_Npart5_randomwalk_PBC_R2_basic_mc' #'Nsteps20000_Nreal100_Nblocks3_intsigma2_intd10_Npart5_pot_exp6.000000_sigma1.000000_factor1.000000_PBC_R2_basic_mc' #'Nsteps20000_Nreal100_Npart5_hardpot_R2_basic_mc'

plotname_totwalk  = infilename_totwalk + '_step_vs_R2' # '.png'
plotname_sections = infilename_sections + '_step_vs_R2'

plotname_totwalk_sparser  = infilename_totwalk + '_step_vs_R2_sparser'
plotname_sections_sparser = infilename_sections + '_step_vs_R2_sparser'

if printall==False:
    plotname_totwalk    = plotname_totwalk    + '_printevery%i.png' % printevery
    plotname_sections   = plotname_sections   + '_printevery%i.png' % printevery
    infilename_totwalk  = infilename_totwalk  + '_printevery%i' % printevery
    infilename_sections = infilename_sections + '_printevery%i' % printevery
    plotname_totwalk_sparser  = plotname_totwalk_sparser  +  '_printevery%i.png' % printevery # Not have sparser plots if we only print sometimes? Or nice to have the option?
    plotname_sections_sparser = plotname_sections_sparser +  '_printevery%i.png' % printevery
else:
    plotname_totwalk    = plotname_totwalk   + '.png'
    plotname_sections   = plotname_sections  + '.png'
    plotname_totwalk_sparser  = plotname_totwalk_sparser  + '.png'
    plotname_sections_sparser = plotname_sections_sparser + '.png'

if nearwall_rw==True or nearwall_mc==True or onedim==True: # Totally different file convention for these ones.
    difftype = ''
    if nearwall_rw==True:
        difftype = 'rw'
    else:
        difftype = 'mc'
    infilename_totwalk  = 'nearwall'+difftype + '_R2_Nsteps%i_Nreal%i_maxstartdist%i_printevery%i' % (Nsteps, Nreal, maxstartdist, printevery)
    infilename_sections = 'nearwall'+difftype + '_R2_Nsteps%i_Nreal%i_maxstartdist%i_Npart%i_printevery%i' % (Nsteps, Nreal, maxstartdist, Nsect, printevery)
    if onedim==True:
        infilename_totwalk  = 'rw1D_R2_Nsteps%i_Nreal%i_printevery%i' % (Nsteps, Nreal, printevery)
        infilename_sections = 'rw1D_R2_Nsteps%i_Nreal%i_Npart%i_printevery%i' % (Nsteps, Nreal, Nsect, printevery)
    plotname_totwalk    = infilename_totwalk  + '.png'
    plotname_sections   = infilename_sections + '.png'
    plotname_totwalk_sparser  = infilename_totwalk + '_step_vs_R2_sparser.png'
    plotname_sections_sparser = infilename_sections + '_step_vs_R2_sparser.png'



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

# Dense plot (this is too dense)
plt.figure(figsize=(6,5))
plt.errorbar(step, R2, yerr=R2_stdv, fmt="none", capsize=2)
plt.xlabel(r'Step', fontsize=16)
plt.ylabel(r'<$R^2$>', fontsize=16)
plt.tight_layout(pad=3.0)#, w_pad=0.0, h_pad=0.5)
#plt.legend(loc="upper right")
plt.title(r'<$R^2$> vs step', fontsize=16) # (last chain)', fontsize=16)
if savefig==True:
    plt.savefig(plotname_totwalk)
else:
    plt.show()


# Sparse plot (experiment with this)
sparsity  = 250
if printall==False:
    sparsity = 3 # Lower if we don't have that many points written to file.
step_sparse    = step[::sparsity]
R2_sparse      = R2[::sparsity]
R2_stdv_sparse = R2_stdv[::sparsity]

plt.figure(figsize=(6,5))
plt.errorbar(step_sparse, R2_sparse, yerr=R2_stdv_sparse, fmt="none", capsize=2)
plt.plot(step_sparse, R2_sparse)
plt.xlabel(r'Step', fontsize=16)
plt.ylabel(r'<$R^2$>', fontsize=16)
plt.tight_layout(pad=3.0)#, w_pad=0.0, h_pad=0.5)
#plt.legend(loc="upper right")
plt.title(r'<$R^2$> vs step, every %i th value' % sparsity, fontsize=16) # (last chain)', fontsize=16)
if savefig==True:
    plt.savefig(plotname_totwalk_sparser)
else:
    plt.show()

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

plt.figure(figsize=(6,5))
for i in range(Nsect):
    plt.errorbar(step_all[i], R2_all[i], yerr=R2_stdv_all[i], fmt="none", capsize=2, label='Section %i' % i)
plt.xlabel(r'Step', fontsize=16)
plt.ylabel(r'<$R^2$>', fontsize=16)
plt.tight_layout(pad=3.0)#, w_pad=0.0, h_pad=0.5)
plt.legend(loc="lower right")
plt.title(r'<$R^2$> vs step', fontsize=16) # (last chain)', fontsize=16)
if savefig==True:
    plt.savefig(plotname_sections)
else:
    plt.show()


plt.figure(figsize=(6,5))
for i in range(Nsect):
    plt.plot(step_all[i], R2_all[i], label='Section %i' % i)
plt.xlabel(r'Step', fontsize=16)
plt.ylabel(r'<$R^2$>', fontsize=16)
plt.tight_layout(pad=3.0)#, w_pad=0.0, h_pad=0.5)
plt.legend(loc="lower right")
plt.title(r'<$R^2$> vs step', fontsize=16) # (last chain)', fontsize=16)
if savefig==False:
    plt.show()


sparsity2 = 10
if printall==False:
    sparsity2 = 2 # Lower if we don't have that many points written to file.
plt.figure(figsize=(6,5))
for i in range(Nsect):
    step_this    = step_all[i]
    R2_this      = R2_all[i]
    R2_stdv_this = R2_stdv_all[i]
    plot_step    = step_this[::sparsity2]
    plot_R2      = R2_this[::sparsity2]
    plot_R2_stdv = R2_stdv_this[::sparsity2]
    plt.errorbar(plot_step, plot_R2, yerr=plot_R2_stdv, fmt="none", capsize=2, label='Section %i' % i)
plt.xlabel(r'Step', fontsize=16)
plt.ylabel(r'<$R^2$>', fontsize=16)
plt.tight_layout(pad=3.0)#, w_pad=0.0, h_pad=0.5)
plt.legend(loc="lower right")
plt.title(r'<$R^2$> vs step, every %i th value' % sparsity2, fontsize=16) # (last chain)', fontsize=16)
if savefig==True:
    plt.savefig(plotname_sections_sparser)
else:
    plt.show()




