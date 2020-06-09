from pylab import *
import matplotlib.pyplot as plt                     # To plot
import numpy as np
import math

saveit = False

damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
spacing = 2
psigma  = 1
print('spacing:', spacing)
print('psigma:', psigma)
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

basepath    = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp+str(psigma) + '/'
inlocation  = basepath + 'log-files/'
outlocation = basepath + 'Test_energy/'

namebase   = 'confignr2_printevery10'
plotname   = outlocation+'d'+str(spacing)+'_'+namebase

infilename = inlocation+'log.'+namebase
infile = open(infilename, 'r')
lines  = infile.readlines()


counter = 0
startit = 0
timestep = []
pot_en   = []
for line in lines:
    words = line.split()
    if startit==1:
        if words[0]=='Loop' or words[0]=='WARNING:':
            break
        timestep.append(int(words[0]))
        pot_en.append(float(words[1]))
    if len(words)>1:
        if words[0]=='Step':
            startit = 1

timestep = np.array(timestep)
pot_en   = np.array(pot_en)

plt.figure(figsize=(6,5))
plt.plot(timestep,pot_en)
plt.xlabel('Time step')
plt.ylabel('Potential energy')
plt.title('Potential energy of bead')
plt.tight_layout()
if saveit==True:
    plt.savefig(plotname)
else:
    plt.show()