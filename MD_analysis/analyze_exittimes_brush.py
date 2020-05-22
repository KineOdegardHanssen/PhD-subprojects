from pylab import *
import matplotlib.pyplot as plt                     # To plot
import numpy as np
import math

damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
spacing = 100
psigma  = 1
print('spacing:', spacing)
print('psigma:', psigma)

savethefig = True

# I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T            = 3
confignrs    = np.arange(1,1001)#300)#20)#101)#1001)#22) 
Nseeds       = len(confignrs)      # So that I don't have to change that much
Nsteps       = 2001 # 20001 before
unitlength   = 1e-9
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime#*writeevery
Npartitions  = 5 # For extracting more walks from one file (but is it really such a random walk here...?)
minlength    = int(floor(Nsteps/Npartitions)) # For sectioning the data
print('timestepsize:', timestepsize)

endlocation          = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing%i/damp%i_diffseedLgv/Brush/Sigma_bead_' % (spacing,damp)+str(psigma) + '/'
filestext            = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])
storelocation        = endlocation + 'Exitanalysis/'

# Plot names (in case we store the plots)
plotname_d2s               = storelocation + 'd2s.png'
plotname_ztrajs_exitvsnot  = storelocation + 'ztrajs_exitvsnot.png'
plotname_ztrajs_average    = storelocation + 'ztrajs_average.png'
plotname_counters          = storelocation + 'counters.png'
plotname_histogram         = storelocation + 'histogram.png'
plotname_histogram_cumul   = storelocation + 'histogram_cumul.png'

# Text files
infilename_ds       = endlocation+'av_ds_'+filestext+'.txt'                        #'lammpsdiffusion_qdrgr_'+namebase+filestext+'_av_ds.txt'
infilename_cuts     = endlocation+'cuts.txt'
infilename_noexits  = endlocation+'noexitzs.txt'

# Plots
plotname             = endlocation+'par_ort_beginning_'+filestext+'_wcuttimes.png'                                 #'lammpsdiffusion_qdrgr_'+namebase+filestext+'.png'

infile_ds = open(infilename_ds, 'r')
firstline = infile_ds.readline()

#Time step; Time; <R^2>; <dx^2>; <dy^2>; <dz^2>; <dx^2+dy^2>
#  [0]       [1]   [2]     [3]    [4]      [5]       [6]
timesteps = []
times     = []
dR2       = []
dx2       = []
dy2       = []
dz2       = []
dpar2     = []

lines = infile_ds.readlines()

for line in lines:
    words = line.split()
    timesteps.append(int(words[0]))
    times.append(float(words[1]))
    dR2.append(float(words[2]))
    dx2.append(float(words[3]))
    dy2.append(float(words[4]))
    dz2.append(float(words[5]))
    dpar2.append(float(words[6]))
infile_ds.close()

timesteps = np.array(timesteps)
times     = np.array(times)
dR2       = np.array(dR2)
dx2       = np.array(dx2)
dy2       = np.array(dy2)
dz2       = np.array(dz2)
dpar2     = np.array(dpar2)

# For generating averages:
# Accumulation
averagez_noexit = np.zeros(Nsteps)
averagez_exit   = np.zeros(Nsteps)
averagez_all    = np.zeros(Nsteps)
# Counter
counter_noexit  = np.zeros(Nsteps)
counter_exit    = np.zeros(Nsteps)

# Data from cut trajectories
infile_cuts = open(infilename_cuts, 'r')
firstline   = infile_cuts.readline()

#Time step; Time; <R^2>; <dx^2>; <dy^2>; <dz^2>; <dx^2+dy^2>
#  [0]       [1]   [2]     [3]    [4]      [5]       [6]
configs  = []
cuttimes = []
runtimes = []
ztrajs   = []
runtimes_flat = []
ztrajs_flat   = []

lines = infile_cuts.readlines()

for line in lines:
    words = line.split()
    lenwords = len(words)
    configs.append(int(words[0]))
    cuttimes.append(float(words[1]))
    # z-coords.:
    ztrajs_this   = []
    runtimes_this = []
    for i in range(2,len(words)):
        zthis = float(words[i])
        thisi = i-2
        ztrajs_this.append(zthis)
        runtimes_this.append(thisi)
        ztrajs_flat.append(zthis)
        runtimes_flat.append(thisi)
        # For averages:
        averagez_exit[thisi] += zthis
        counter_exit[thisi]  += 1
    runtimes.append(runtimes_this)
    ztrajs.append(ztrajs_this)
infile_cuts.close()

configs  = np.array(configs)  # Do I even need the configs?
cuttimes = np.array(cuttimes)

## Data from uncut trajectories
infile_noexits = open(infilename_noexits, 'r')
firstline      = infile_noexits.readline()

configs_noexit = []
runtimes_noexit = []
ztrajs_noexit   = []
runtimes_flat_noexit = []
ztrajs_flat_noexit   = []

lines = infile_noexits.readlines()

for line in lines:
    words = line.split()
    lenwords = len(words)
    configs_noexit.append(int(words[0]))
    # z-coords.:
    ztrajs_this   = []
    runtimes_this = []
    for i in range(1,len(words)):
        zthis = float(words[i])
        thisi = i-1
        ztrajs_this.append(zthis)
        runtimes_this.append(thisi)
        ztrajs_flat_noexit.append(zthis)
        runtimes_flat_noexit.append(thisi)
        # For averages:
        averagez_noexit[thisi] += zthis
        counter_noexit[thisi]  += 1      # Not sure I really need this, but...
    runtimes_noexit.append(runtimes_this)
    ztrajs_noexit.append(ztrajs_this)
infile_cuts.close()

configs_noexit  = np.array(configs_noexit)  # Do I even need the configs?

# Generating averages


for i in range(Nsteps):
    averagez_all[i] = (averagez_exit[i]+averagez_noexit[i])/(counter_noexit[i]+counter_exit[i])
    if counter_exit[i]!=0:
        averagez_exit[i]   /= counter_exit[i]
    if counter_noexit[i]!=0:
        averagez_noexit[i] /= counter_noexit[i]

## To determine the range
plt.figure(figsize=(6,5))
plt.plot(timesteps, dR2, label=r'$<R^2>$')
plt.plot(timesteps, dx2, label=r'$<dx^2>$')
plt.plot(timesteps, dy2, label=r'$<dy^2>$')
plt.plot(timesteps, dz2, label=r'$<dz^2>$')
plt.plot(timesteps, dpar2, label=r'$<dx^2+dy^2>$')
plt.xlabel(r'Index (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
if savethefig==True:
    plt.savefig(plotname_d2s)


plt.figure(figsize=(6,5))
plt.plot(runtimes_flat, ztrajs_flat, ',', label=r'Exits')
plt.plot(runtimes_flat_noexit, ztrajs_flat_noexit, ',', label=r'Does not exit')
plt.xlabel(r'Index (s)')
plt.ylabel(r'$dz^2$ [in unit length]')
plt.title(r'$dz^2$, exit vs no exit, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
if savethefig==True:
    plt.savefig(plotname_ztrajs_exitvsnot)

plt.figure(figsize=(6,5))
plt.plot(timesteps, averagez_exit, label=r'Exits')
plt.plot(timesteps, averagez_noexit, label=r'Does not exit')
plt.plot(timesteps, averagez_all, label=r'All')
plt.xlabel(r'Index (s)')
plt.ylabel(r'$<dz^2>$ [in unit length]')
plt.title(r'$<dz^2>$, exit vs no exit, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
if savethefig==True:
    plt.savefig(plotname_ztrajs_average)

plt.figure(figsize=(6,5))
plt.plot(timesteps, counter_exit, label=r'Exits')
#plt.plot(timesteps, counter_noexit, label=r'Does not exit')
plt.xlabel(r'Index (s)')
plt.ylabel(r'Number of configs')
plt.title(r'Configs exited, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
if savethefig==True:
    plt.savefig(plotname_counters)

# Saving
'''
xmax_plot = 100
plt.figure(figsize=(6,5))
plt.plot(timesteps, dR2, label=r'$<R^2>$')
plt.plot(timesteps, dx2, label=r'$<dx^2>$')
plt.plot(timesteps, dy2, label=r'$<dy^2>$')
plt.plot(timesteps, dz2, label=r'$<dz^2>$')
plt.plot(timesteps, dpar2, label=r'$<dx^2+dy^2>$')
plt.xlabel(r'Index (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.axis([0, xmax_plot, 0, max(dR2[0:xmax_plot])])
plt.savefig(plotname)
'''

plt.figure(figsize=(6,5))
plt.hist(cuttimes, bins=100)
plt.xlabel(r'Exit time (Time of cut)')
plt.ylabel(r'Number of exits')
plt.title(r'Histogram of exit times, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
if savethefig==True:
    plt.savefig(plotname_histogram)

plt.figure(figsize=(6,5))
plt.hist(cuttimes, bins=100, cumulative=True)
plt.xlabel(r'Exit time (Time of cut)')
plt.ylabel(r'Number of exits')
plt.title(r'Cumulative histogram of exit times, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
if savethefig==True:
    plt.savefig(plotname_histogram_cumul)

if savethefig==False:
    plt.show()



