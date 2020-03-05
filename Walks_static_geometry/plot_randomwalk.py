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

Nsteps = 20000
Nsect  = 100

### Setting file names
infilename_totwalk  = 'Nsteps20000_Nreal100_hardpot_R2_basic_mc'
infilename_sections = 'Nsteps20000_Nreal100_Npart5_hardpot_R2_basic_mc'

plotname_totwalk  = infilename_totwalk + '_step_vs_R2.png'
plotname_sections = infilename_sections + '_step_vs_R2.png'

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
plt.savefig(plotname_totwalk)

# Sparse plot (experiment with this)
sparsity  = 250
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
plt.savefig(plotname_totwalk_sparser)

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
plt.savefig(plotname_sections)

sparsity2 = 50
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
plt.savefig(plotname_sections_sparser)



