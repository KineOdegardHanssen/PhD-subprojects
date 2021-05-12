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

plottest = False # True # 

Nrms = 10 # Number of averages 

DR_estimates = np.zeros(Nrms)
Dz_estimates = np.zeros(Nrms)
Dpar_estimates = np.zeros(Nrms)
# Constant terms (for plotting purposes)
bR_estimates = np.zeros(Nrms)
bz_estimates = np.zeros(Nrms)
bpar_estimates = np.zeros(Nrms)

def rmsd(x,y):
    Nx = len(x)
    Ny = len(y)
    if Nx!=Ny:
        print('WARNING! Nx!=Ny. Could not calculate rmsd value')
        return 'WARNING! Nx!=Ny. Could not calculate rmsd value'
    delta = 0
    for i in range(Nx):
        delta += (x[i]-y[i])*(x[i]-y[i])
    delta = np.sqrt(delta/(Nx-1))
    return delta

def avg_and_rms(x):
    N = len(x)
    avgx = np.mean(x)
    print('avgx:',avgx)
    rmsx = 0
    for i in range(N):
        rmsx += (avgx-x[i])**2
    rmsx = np.sqrt(rmsx/(N-1))
    return avgx,rmsx

# Settings
Nconfigs    = 1000
Nplacements = 11
#
damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
popup_plots = False
testmode = False
long    = False #True
spacing = 7
psigma  = 1
density = 0.238732414637843 # Yields mass 1 for bead of radius 1 nm
print('spacing:', spacing)
print('psigma:', psigma)
# I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T        = 3
plotseed = 0
plotdirs = False
test_sectioned = False
zhigh          = 250
zlow           = -50
confignrs      = np.arange(1,Nconfigs+1)
beadplacements = np.arange(Nplacements-1,Nplacements+1) # DEPRECATED
Nconfigs       = len(confignrs)      # So that I don't have to change that much
Nplacements    = len(beadplacements)
Nfiles         = Nconfigs*Nplacements
N_per_rms      = int(Nconfigs/Nrms)
filescounter   = 0
if long==True:
    Nsteps     = 10001
else:
    Nsteps     = 2001
#writeevery    = 10                  # I'm writing to file every this many time steps
unitlength     = 1e-9
unittime       = 2.38e-11 # s
timestepsize   = 0.00045*unittime#*writeevery
Npartitions    = 5 # For extracting more walks from one file (but is it really such a random walk here...?)
minlength      = int(floor(Nsteps/Npartitions)) # For sectioning the data
gatheredconfignrs = np.arange(1,Nfiles+1) # For some old printing
print('timestepsize:', timestepsize)

basepath        = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/Spacing'+str(spacing)+'/'
location_config = basepath + 'Radius1/Initial_configs/Before_bead/'
basepath        = basepath + 'Radius' + str(psigma) + '/'
endlocation     = basepath + 'Nocut/'

print('End location:', endlocation)

filestext   = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])+'_placements'+str(beadplacements[0])+'to'+str(beadplacements[-1])

# Text files
outfilename = endlocation+'diffusion'+filestext+'_better_rms_Nestimates%i.txt' % Nrms
# Plots
plotname    = endlocation+'_nocut_better_rms_Nestimates%i_compare.png' % Nrms


## Getting interval to perform the fit:
rangefilename = endlocation+'indices_for_fit_better_rms.txt'
rangefile     = open(rangefilename,'r')
line          = rangefile.readline()   
startindex    = int(line.split()[0])
endindex      = int(line.split()[1])
rangefile.close()


# These are all squared:
# All together:
allRs     = []
alldxs    = []
alldys    = []
alldzs    = []
# Averages:
tot_averageR2s  = np.zeros(Nsteps)   # Distances, squared
tot_averagedx2s = np.zeros(Nsteps)
tot_averagedy2s = np.zeros(Nsteps)
tot_averagedz2s = np.zeros(Nsteps)
tot_averagedparallel2 = np.zeros(Nsteps) # Distances, squared
tot_average_counter   = np.zeros(Nsteps)

Nins          = []

linestart_data = 22
skippedfiles = 0
unphysical   = 0

for j in range(Nrms):
    configbase = N_per_rms*j+1
    averageR2s  = np.zeros(Nsteps)   # Distances, squared
    averagedx2s = np.zeros(Nsteps)
    averagedy2s = np.zeros(Nsteps)
    averagedz2s = np.zeros(Nsteps)
    averagedxs  = np.zeros(Nsteps)   # Distances, non-squared
    averagedys  = np.zeros(Nsteps)
    averagedzs  = np.zeros(Nsteps)
    averagedparallel2 = np.zeros(Nsteps) # Distances, squared
    average_counter   = np.zeros(Nsteps)
    for k in range(N_per_rms):
        confignr = k+configbase
        print('On config number:', confignr)
   
        if long==True:
            infilename_free = basepath+'long/'+'freeatom_confignr'+str(confignr)+'_beadplacement10_long.lammpstrj'
        else:
            infilename_free = basepath+'freeatom_confignr'+str(confignr)+'_beadplacement10.lammpstrj'
        
        try:
            infile_free = open(infilename_free, "r")
        except:
            try:
                if long==True:
                    infilename_free = basepath+'long/'+'freeatom_confignr'+str(confignr)+'_beadplacement11_long.lammpstrj'
                else:
                   infilename_free = basepath+'freeatom_confignr'+str(confignr)+'_beadplacement11.lammpstrj'
                infile_free = open(infilename_free, "r")
            except:
                print('file name that failed:', infilename_free)
                skippedfiles += 1
                continue # Skipping this file if it does not exist
        filescounter += 1
        
        # Moving on, if the file exists
        lines = infile_free.readlines() # This takes some time
        # Getting the number of lines, etc.
        totlines = len(lines)         # Total number of lines
        lineend = totlines-1          # Index of last element
        
        # Extracting the number of atoms:
        words = lines[3].split()
        Nall = int(words[0])
        N    = Nall
        
        skiplines   = 9             # If we hit 'ITEM:', skip this many steps...
        skipelem    = 0
        sampleevery = 0
        i           = int(math.ceil(skipelem*(Nall+9)))
        skiplines  += (Nall+skiplines)*sampleevery # Check!
        
        # Setting arrays for treatment:
        
        # Setting arrays for treatment:
        vx = [] # x-component of the velocity. 
        vy = []
        vz = []
        positions = [] # Fill with xu,yu,zu-coordinates. One position for each time step. # We have the same time step for the whole of the simulation. Can bother with time later
        times     = [] # Not useful anymore?
        
        # For now: Only find the largest z-position among the beads in the chain. # Otherwise I must read off the atom number (will do that with the free bead anyways.)
        maxz = -1000 # No z-position is this small.
        
        time_start = time.process_time()
        counter = 0
        zs_fortesting = []
        while i<totlines:
            words = lines[i].split()
            if words[0]=='ITEM:':
                if words[1]=='TIMESTEP':
                    words2 = lines[i+1].split() # The time step is on the next line
                    t = float(words2[0])
                    times.append(t)
                    i+=skiplines
                elif words[1]=='NUMBER': # These will never kick in. 
                    i+=7
                elif words[1]=='BOX':
                    i+=5
                elif words[1]=='ATOMS':
                    i+=1
            elif len(words)<8:
                i+=1
            else:
                # Find properties
                # Order:  id  type mol ux  uy  uz  vx  vy   vz
                #         [0] [1]  [2] [3] [4] [5] [6] [7]  [8]
                ind      = int(words[0])-1 # Atom ids go from zero to N-1.
                #atomtype = int(words[1]) 
                #molID    = int(words[2])
                x        = float(words[3])
                y        = float(words[4])
                z        = float(words[5])
                positions.append(np.array([x,y,z]))
                vx.append(float(words[6]))
                vy.append(float(words[7]))
                vz.append(float(words[8]))
                zs_fortesting.append(z)
                counter+=1
                i+=1
        infile_free.close()
        dt = (times[1]-times[0])*timestepsize # This might be handy
        #print('dt:', dt)
        #print('times[1]-times[0]:',times[1]-times[0])
        times_single = np.arange(Nsteps)#*dt
        times_single_real = np.arange(Nsteps)*dt
        
        time_end = time.process_time()
        
        if max(zs_fortesting)>zhigh:
            if testmode==True:
                print('Cut, max')
            unphysical+=1
            continue
        if min(zs_fortesting)<zlow:
            if testmode==True:
                print('Cut, min')
            unphysical+=1
            continue
        
        if len(zs_fortesting)!=Nsteps: # Something is wrong with the file. Skip
            print('Warning! File error!')
            skippedfiles+=1
            continue
        Nin   = Nsteps
        
        #######
        # Will divide into several walks with different starting points later on
        # Should store for RMS too? ... Then I need to have more starting points or more time frames.
        
        # Finding R2s and corresponding times
        # Maybe... Not plot the first ones since it takes some time to equilibrate.
        # But how much do I need to cut? Varies from plot to plot how much makes sense to cut (if it is even possible to say...)
        # Go for some percentage of all points?
        R_temp  = []
        dx_temp = []
        dy_temp = []
        dz_temp = []
        step_temp = []
        
        R_temp.append(0)
        dx_temp.append(0)
        dy_temp.append(0)
        dz_temp.append(0)
        step_temp.append(0)
        startpos = positions[0]
        for i in range(1,Nin):       
            this = positions[i]
            dist = this-startpos
            dx   = dist[0]         # Signed
            dy   = dist[1]
            dz   = dist[2]
            R2   = np.dot(dist,dist) # Squared
            dx2  = dx*dx
            dy2  = dy*dy
            dz2  = dz*dz
            # All together:
            allRs.append(R2)
            alldxs.append(dx2)
            alldys.append(dy2)
            alldzs.append(dz2)
            # Averages:
            averageR2s[i] +=R2   # Distance
            averagedx2s[i]+=dx2
            averagedy2s[i]+=dy2
            averagedz2s[i]+=dz2
            averagedxs[i] +=dx   # Distance, signed
            averagedys[i] +=dy
            averagedzs[i] +=dz
            averagedparallel2[i] += dx2+dy2 # Distance
            average_counter[i] +=1
            # Total averages: ############### I only need these for plotting
            tot_averageR2s[i] +=R2   # Distances, squared
            tot_averagedx2s[i] +=dx2
            tot_averagedy2s[i] +=dy2
            tot_averagedz2s[i] +=dz2
            tot_averagedparallel2[i] += dx2+dy2 # Distances, squared
            tot_average_counter[i]+=1
            ################################
        
        
        Nins.append(Nin)


    Ninbrush = Nsteps # Default, in case it does not exit the brush
    for i in range(1,Nsteps):
        counter = average_counter[i]
        if counter!=0:
            # Distance squared
            averageR2s[i]/=counter
            averagedx2s[i]/=counter
            averagedy2s[i]/=counter
            averagedz2s[i]/=counter
            averagedparallel2[i]/= counter
            # Distance, signed:
            averagedxs[i]/=counter
            averagedys[i]/=counter
            averagedzs[i]/=counter
        else:
            Ninbrush = i-1
            break
    
    times_single      = times_single[0:Ninbrush]
    times_single_real = times_single_real[0:Ninbrush]
    # Distance squared
    averageR2s  = averageR2s[0:Ninbrush]
    averagedx2s = averagedx2s[0:Ninbrush]
    averagedy2s = averagedy2s[0:Ninbrush]
    averagedz2s = averagedz2s[0:Ninbrush]
    averagedparallel2 = averagedparallel2[0:Ninbrush]
    # Distance to the first power
    averagedxs = averagedxs[0:Ninbrush]
    averagedys = averagedys[0:Ninbrush]
    averagedzs = averagedzs[0:Ninbrush]

    print('len(averageR2s):',len(averageR2s))
    print('Nsteps:', Nsteps)
    
    ## Average, SI units:
    # Distance:
    averageR2s_SI  = averageR2s*unitlength**2
    averagedx2s_SI = averagedx2s*unitlength**2
    averagedy2s_SI = averagedy2s*unitlength**2
    averagedz2s_SI = averagedz2s*unitlength**2
    averagedparallel2_SI = averagedparallel2*unitlength**2
    
    #print('----------------------------------')
    #print(' averageR2s_SI[2]:', averageR2s_SI[2])
    
    # Insert interval
    coeffs_poly, covs = polyfit(times_single_real[startindex:endindex], averageR2s_SI[startindex:endindex], 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
    a_poly_SI = coeffs_poly[0]
    b_poly_SI = coeffs_poly[1]
    rms_D_poly_SI = np.sqrt(covs[0,0])/6.
    rms_b_poly_SI = np.sqrt(covs[1,1])
    D_poly_SI = a_poly_SI/6.
    
    DR_estimates[j] = D_poly_SI
    bR_estimates[j] = b_poly_SI
    
    fit_poly_SI = a_poly_SI*times_single_real+b_poly_SI
    if plottest==True:
       plt.figure(figsize=(6,5))
       plt.plot(times_single_real, averageR2s_SI)
       plt.plot(times_single_real, fit_poly_SI, '--')
       plt.xlabel('Time [s]')
       plt.ylabel('$<dR^2>$')
       plt.title('Test, R')
       plt.show()
    
    coeffs_poly, covs = polyfit(times_single_real[startindex:endindex], averagedparallel2_SI[startindex:endindex], 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
    a_poly_SI = coeffs_poly[0]
    b_poly_SI = coeffs_poly[1]
    rms_D_poly_SI = np.sqrt(covs[0,0])/6.
    rms_b_poly_SI = np.sqrt(covs[1,1])
    D_poly_SI = a_poly_SI/4.
    
    Dpar_estimates[j] = D_poly_SI
    bpar_estimates[j] = b_poly_SI
    
    fit_poly_SI = a_poly_SI*times_single_real+b_poly_SI
    if plottest==True:
       plt.figure(figsize=(6,5))
       plt.plot(times_single_real, averagedparallel2_SI)
       plt.plot(times_single_real, fit_poly_SI, '--')
       plt.xlabel('Time [s]')
       plt.ylabel('$<dparallel^2>$')
       plt.title('Test, parallel')
       plt.show()
    
    coeffs_poly, covs = polyfit(times_single_real[startindex:endindex], averagedz2s_SI[startindex:endindex], 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
    a_poly_SI = coeffs_poly[0]
    b_poly_SI = coeffs_poly[1]
    rms_D_poly_SI = np.sqrt(covs[0,0])/6.
    rms_b_poly_SI = np.sqrt(covs[1,1])
    D_poly_SI = a_poly_SI/2.
    
    Dz_estimates[j] = D_poly_SI
    bz_estimates[j] = b_poly_SI
    
    fit_poly_SI = a_poly_SI*times_single_real+b_poly_SI
    if plottest==True:
       plt.figure(figsize=(6,5))
       plt.plot(times_single_real, averagedz2s_SI)
       plt.plot(times_single_real, fit_poly_SI, '--')
       plt.xlabel('Time [s]')
       plt.ylabel('$<dz^2>$')
       plt.title('Test, z')
       plt.show()
    
for i in range(1,Nsteps):
    counter = tot_average_counter[i]
    if counter!=0:
        # Distance squared
        tot_averageR2s[i]/=counter
        tot_averagedx2s[i]/=counter
        tot_averagedy2s[i]/=counter
        tot_averagedz2s[i]/=counter
        tot_averagedparallel2[i]/= counter
    else:
        Ninbrush = i-1
        break

tot_averageR2s_SI  = tot_averageR2s*unitlength**2
tot_averagedx2s_SI = tot_averagedx2s*unitlength**2
tot_averagedy2s_SI = tot_averagedy2s*unitlength**2
tot_averagedz2s_SI = tot_averagedz2s*unitlength**2
tot_averagedparallel2_SI = tot_averagedparallel2*unitlength**2


DR_avg, DR_rms = avg_and_rms(DR_estimates)
Dz_avg, Dz_rms = avg_and_rms(Dz_estimates)
Dpar_avg, Dpar_rms = avg_and_rms(Dpar_estimates)
# Constant term for plotting purposes:
bR_avg, bR_rms = avg_and_rms(bR_estimates)
bz_avg, bz_rms = avg_and_rms(bz_estimates)
bpar_avg, bpar_rms = avg_and_rms(bpar_estimates)

fit_DR_SI   = 6*DR_avg*times_single_real+bR_avg
fit_Dz_SI   = 2*Dz_avg*times_single_real+bz_avg
fit_Dpar_SI = 4*Dpar_avg*times_single_real+bpar_avg
   
outfile = open(outfilename,'w')
outfile.write('D_R2  sigmaD_R2; D_z2  sigmaD_z2; D_par2 sigmaD_par2\n')
outfile.write('%.5e %.5e %.5e %.5e %.5e %.5e\n' % (DR_avg, DR_rms, Dz_avg, Dz_rms, Dpar_avg, Dpar_rms))
outfile.close()

print('len(positions):',len(positions))
print('len(times_single):',len(times_single))

# Plot 
plt.figure(figsize=(6,5))
plt.plot(times_single_real, tot_averageR2s_SI, label=r'$<dR^2>$')
plt.plot(times_single_real, tot_averagedz2s_SI, label=r'$<dz^2>$')
plt.plot(times_single_real, tot_averagedparallel2_SI, label=r'$<dx^2+dy^2>$')
plt.plot(times_single_real, fit_DR_SI, '--', label=r'$<dR^2>$, fit')
plt.plot(times_single_real, fit_Dz_SI, '--', label=r'$<dz^2>$, fit')
plt.plot(times_single_real, fit_Dpar_SI, '--', label=r'$<dx^2+dy^2>$, fit')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'Averaged RMSD in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname)


print('counters:',tot_average_counter)

print('spacing:', spacing)
print('psigma:', psigma)

print('Testing/error search:')
print('average_counter[1]:', average_counter[1])
'''
for i in range(len(average_counter)):
    print(i,' average_counter:', average_counter[i])
'''
print('unphysical:',unphysical)
plt.show()