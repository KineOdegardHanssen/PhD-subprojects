import matplotlib.pyplot as plt
import numpy as np


# Resolution
lx = 0.5  # nm. Want this smaller than sigma.
ly = lx
thr = 0.002 # Revisit this after first run

damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
spacings = [1,2,5,8,100]
psigma  = 1
print('spacings:', spacings)
print('psigma:', psigma)
# I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T        = 3
plotseed = 0
plotdirs = False
test_sectioned = False
confignrs    = np.arange(1,1001)
Nseeds       = len(confignrs)      # So that I don't have to change that much
maxz_av      = 0
filescounter = 0
Nsteps       = 2001 # 20001 before
#writeevery   = 10                  # I'm writing to file every this many time steps
unitlength   = 1e-9
unittime     = 2.38e-11 # s
timestepsize = 0.00045*unittime#*writeevery
print('timestepsize:', timestepsize)

filestext            = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])
outlocation          = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_forest/Projected_trajectories/'
plotfilename         = outlocation +'midlines_lx'+str(lx)+'_nocut.png'

plt.figure(figsize=(6,5))
for spacing in spacings:
    xpos = []
    intensity = []
    endlocation_in       = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_forest/Spacing'+str(spacing)+'/'
    endlocation = endlocation_in + 'Results/' # Nocut is default
    # Should I have printed these somewhere else?
    infilename  = endlocation+'midline_'+filestext+'_d'+str(spacing)+'_lx'+str(lx)+'_nocut.txt'
    infile      = open(infilename,'r')
    lines = infile.readlines()
    for line in lines:
        words = line.split()
        if len(words)!=0:
            xpos.append(float(words[0]))
            intensity.append(float(words[1])/2001) # Occurences per time step per file
    infile.close()
    plt.plot(xpos,intensity,label='d=%i' %spacing)

plt.xlabel(r'$x/d$')
plt.ylabel(r'Probability') # Is this correct?!
plt.legend(loc='upper right') # Do I need an outside box for this?
plt.title('x-position vs probability')
plt.tight_layout()
plt.savefig(plotfilename)
plt.show()

