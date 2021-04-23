import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob

def avg_and_rms(x):
    N = len(x)
    avgx = np.mean(x)
    print('avgx:',avgx)
    rmsx = 0
    for i in range(N):
        rmsx += (avgx-x[i])**2
    rmsx = np.sqrt(rmsx/(N-1))
    return avgx,rmsx

ljfactor   = 10
Nintervals = 10 # F.ex.

# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
bulkdiffusion = False #False
substrate     = False

#spacing = 7
spacings = [2.5]#[25,50,75,100]#[7,8,10,15]#[3,4,5,6]#[1,1.25,1.5,2]#
psigma   = 1         # So far, we need to treat one and one sigma.
damp     = 10

Nsteps = 2001
confignrs = np.arange(1,101)
placements = np.arange(1,11)

# Make file names

for spacing in spacings:
    filestext  = '_config'+str(confignrs[0])+'to'+str(confignrs[-1])+'_placements'+str(placements[0])+'to'+str(placements[-1])
    
    endlocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/Spacing'+str(spacing)+'/Radius' + str(psigma) + '/Nocut/'
    
    infilename      = endlocation+'av_ds'+filestext+'_nocut_ljfactor'+str(ljfactor)+'.txt'
    outfilename     = endlocation+'diffusion'+filestext+'_nocut_better_rms_ljfactor'+str(ljfactor)+'_Nintervals%i.txt' % Nintervals # '_nocut' here too?
    metaname        = endlocation+'diffusion_metadata'+filestext+'_nocut_better_rms_ljfactor'+str(ljfactor)+'_Nintervals%i.txt' % Nintervals
    plotname_R      = endlocation+'diffusion_R'+filestext+'_nocut_better_rms_ljfactor'+str(ljfactor)+'_Nintervals%i.png' % Nintervals
    plotname_dz     = endlocation+'diffusion_dz'+filestext+'_nocut_better_rms_ljfactor'+str(ljfactor)+'_Nintervals%i.png' % Nintervals
    plotname_dpar   = endlocation+'diffusion_dpar'+filestext+'_nocut_better_rms_ljfactor'+str(ljfactor)+'_Nintervals%i.png' % Nintervals
    
    # Choosing which part to fit: in file
    rangefilename = endlocation+'indices_for_fit_better_rms_ljfactor'+str(ljfactor)+'.txt'
    rangefile     = open(rangefilename,'r')
    line          = rangefile.readline()
    
    startindex    = int(line.split()[0])
    endindex      = int(line.split()[1])
    rangefile.close()

    # I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

    # Should divide into folders in a more thorough manner?
    # Extracting the correct names (all of them)
    plotseed = 0
    plotdirs = False
    test_sectioned = False
    unitlength   = 1e-9
    unittime     = 2.38e-11 # s
    timestepsize = 0.00045*unittime
    print('timestepsize:', timestepsize)
    Npartitions = 5 # For extracting more walks from one file (but is it really such a random walk here...?)


    # Set up arrays
    Nsteps = endindex-startindex
    times  = np.zeros(Nsteps)
    dR2s   = np.zeros(Nsteps) 
    dz2s   = np.zeros(Nsteps) 
    dpar2s = np.zeros(Nsteps)
    
    DRs   = []
    Dzs   = []
    Dpars = []
    
    time_all    = []
    fit_R_all   = []
    fit_z_all   = []
    fit_par_all = []
    
    # Read file
    infile = open(infilename,'r')
    lines = infile.readlines()

    #    	0			1		2		3		4			5		6
    # (times_single[i], times_single_real[i], averageRs_SI[i], averagedxs_SI[i], averagedys_SI[i], averagedzs_SI[i], averagedparallel_SI[i]))
    for i in range(startindex,endindex):
        j     = i-startindex
        line  = lines[i]
        words = line.split()
        if i>=startindex and i<=endindex: # The hacky way to do it. Could have had three loops...
            times[j]  = float(words[1])
            dR2s[j]   = float(words[2])
            dz2s[j]   = float(words[5])
            dpar2s[j] = float(words[6])

    infile.close()

    # line fit, all
    coeffs, covs = polyfit(times, dR2s, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
    a_R2 = coeffs[0]
    b_R2 = coeffs[1]
    D_R2 = a_R2/6.
    rms_D_R2 = np.sqrt(covs[0,0])/6.
    rms_b_R2 = np.sqrt(covs[1,1])
     
    fit_R2 = a_R2*times+b_R2

    # line fit, dz
    coeffs, covs = polyfit(times, dz2s, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
    a_z2 = coeffs[0]
    b_z2 = coeffs[1]
    D_z2 = a_z2/2.
    rms_D_z2 = np.sqrt(covs[0,0])/2.
    rms_b_z2 = np.sqrt(covs[1,1])

    fit_z2 = a_z2*times+b_z2

    # line fit, parallel
    coeffs, covs = polyfit(times, dpar2s, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
    a_par2     = coeffs[0]
    b_par2     = coeffs[1]
    D_par2     = a_par2/4.
    rms_D_par2 = np.sqrt(covs[0,0])/4.
    rms_b_par2 = np.sqrt(covs[1,1])

    fit_par2 = a_par2*times+b_par2
    
    DRs.append(D_R2)
    Dzs.append(D_z2)
    Dpars.append(D_par2)
    
    
    time_all.append(times)
    fit_R_all.append(fit_R2)
    fit_z_all.append(fit_z2)
    fit_par_all.append(fit_par2)
    
    si = 0
    di = int(math.floor((endindex-startindex)/Nintervals))
    ei = si+di
    for i in range(Nintervals):
        # Do stuff
        # line fit, all
        coeffs, covs = polyfit(times[si:ei], dR2s[si:ei], 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
        a_R2 = coeffs[0]
        b_R2 = coeffs[1]
        D_R2 = a_R2/6.
        rms_D_R2 = np.sqrt(covs[0,0])/6.
        rms_b_R2 = np.sqrt(covs[1,1])
     
        fit_R2 = a_R2*times[si:ei]+b_R2

        # line fit, dz
        coeffs, covs = polyfit(times[si:ei], dz2s[si:ei], 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
        a_z2 = coeffs[0]
        b_z2 = coeffs[1]
        D_z2 = a_z2/2.
        rms_D_z2 = np.sqrt(covs[0,0])/2.
        rms_b_z2 = np.sqrt(covs[1,1])
        
        fit_z2 = a_z2*times[si:ei]+b_z2
        
        # line fit, parallel
        coeffs, covs = polyfit(times[si:ei], dpar2s[si:ei], 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
        a_par2     = coeffs[0]
        b_par2     = coeffs[1]
        D_par2     = a_par2/4.
        rms_D_par2 = np.sqrt(covs[0,0])/4.
        rms_b_par2 = np.sqrt(covs[1,1])
        
        fit_par2 = a_par2*times[si:ei]+b_par2
        
        DRs.append(D_R2)
        Dzs.append(D_z2)
        Dpars.append(D_par2)
        
        time_all.append(times[si:ei])
        fit_R_all.append(fit_R2)
        fit_z_all.append(fit_z2)
        fit_par_all.append(fit_par2)

        # Update indices
        si +=di # Start index
        ei +=di
    
    D_R2,rms_D_R2 = avg_and_rms(DRs)
    D_z2,rms_D_z2 = avg_and_rms(Dzs)
    D_par2,rms_D_par2 = avg_and_rms(Dpars)
    
    outfile = open(outfilename,'w')
    outfile.write('D_R2  sigmaD_R2; D_z2  sigmaD_z2; D_par2 sigmaD_par2\n')
    outfile.write('%.5e %.5e %.5e %.5e %.5e %.5e\n' % (D_R2, rms_D_R2, D_z2, rms_D_z2, D_par2, rms_D_par2))
    outfile.close()

    metafile = open(metaname, 'w')
    metafile.write('startindex: %i\n' % startindex)
    metafile.write('endindex: %i\n' % endindex)
    metafile.close()

    # Plotting for verification
    plt.figure(figsize=(6,5))
    plt.plot(times, dR2s, label='Average, brush')
    for i in range(Nintervals+1):
        plt.plot(time_all[i], fit_R_all[i], '--')
    plt.xlabel(r'Time (s)')
    plt.ylabel(r'Distance$^2$ [in unit length]')
    plt.title('RMSD in bulk, d = %i nm, SI' % spacing)
    plt.tight_layout()
    plt.legend(loc='upper left')
    plt.savefig(plotname_R)

    plt.figure(figsize=(6,5))
    plt.plot(times, dz2s, label='Average, brush')
    for i in range(Nintervals+1):
        plt.plot(time_all[i], fit_z_all[i], '--')
    plt.xlabel(r'Time (s)')
    plt.ylabel(r'$dz^2$ [in unit length]')
    plt.title('$dz^2$ in bulk, d = %i nm, SI' % spacing)
    plt.tight_layout()
    plt.legend(loc='upper left')
    plt.savefig(plotname_dz)

    plt.figure(figsize=(6,5))
    plt.plot(times, dpar2s, label='Average, brush')
    for i in range(Nintervals+1):
        plt.plot(time_all[i], fit_par_all[i], '--')
    plt.xlabel(r'Time (s)')
    plt.ylabel(r'$dx^2+dy^2$ [in unit length]')
    plt.title('$dx^2+dy^2$ in bulk, d = %i nm, SI' % spacing)
    plt.tight_layout()
    plt.legend(loc='upper left')
    plt.savefig(plotname_dpar)

    print('spacing:', spacing)
    print('psigma:', psigma)
    print('damp:', damp)