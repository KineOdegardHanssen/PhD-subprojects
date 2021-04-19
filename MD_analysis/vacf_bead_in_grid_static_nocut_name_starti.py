from skimage import morphology, measure, io, util   # For performing morphology operations (should maybe import more from skimage...)
from mpl_toolkits.mplot3d import Axes3D             # Plotting in 3D
import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
from scipy.ndimage import measurements, convolve    # Should have used the command below, but changing now will lead to confusion
from scipy import ndimage                           # For Euclidean distance measurement
import data_treatment as datr
import numpy as np
import random
import math
import time
import os
import glob
import copy

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

# spacing 1 not here, spacings 3.5, 4.5 long.
spacings = [75,100]#[1.25,1.5,2,2.5,3,4,5,6,7,8,10,15,25,50,75,100]

damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
popup_plots = False
long    = False
psigma  = 1
density = 0.238732414637843 # Yields mass 1 for bead of radius 1 nm
print('psigma:', psigma)

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T        = 3
plotdirs = False
zhigh          = 250
zlow           = -50
confignrs      = np.arange(1,101)
beadplacements = np.arange(1,11)
Nseeds         = len(confignrs)      # So that I don't have to change that much
maxz_av        = 0
filescounter   = 0
if long==True:
    Nsteps     = 10001
else:
    Nsteps     = 2001 # 20001 before
#writeevery    = 10                  # I'm writing to file every this many time steps
unitlength     = 1e-9
unittime       = 2.38e-11 # s
timestepsize   = 0.00045*unittime#*writeevery
filestext      = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])+'_bpl'+str(beadplacements[0])+'to'+str(beadplacements[-1])
print('timestepsize:', timestepsize)
linestart_data = 22
istart = 15 # 2 # For 'equilibration'
starti = istart # Due to stupid naming error

for spacing in spacings:
    endlocation_in       = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/Spacing' + str(spacing)+'/Radius'+str(psigma) + '/'
    endlocation          = endlocation_in +'Nocut/'
    configfolder         = endlocation_in + 'Initial_configs/Before_bead/'
    
    # Text files
    outfilename_vacf     = endlocation+'vacf_'+filestext+'_nocut_starti%i' % starti
    outfilename_Dvacf    = endlocation+'Dvacf_'+filestext+'_nocut_starti%i' % starti
    
    # Plots
    plotname             = endlocation+filestext+'vacf_nocut_starti%i' % starti
    
    if long==True:
        outfilename_vacf     = outfilename_vacf+'_long.txt'
        outfilename_Dvacf    = outfilename_Dvacf+'_long.txt'
        
        # Plots
        plotname             = plotname+'_long.png'
    else:
        outfilename_vacf     = outfilename_vacf+'.txt'
        outfilename_Dvacf    = outfilename_Dvacf+'.txt'
        
        # Plots
        plotname             = plotname+'.png'
        
    ## Setting arrays
    # These are all squared:
    # All together:
    allRs     = []
    average_counter = np.zeros(Nsteps-istart)
    # Velocity autocorrelation function:
    vacf_tot = np.zeros(Nsteps-istart)
    vacf_x   = np.zeros(Nsteps-istart)
    vacf_y   = np.zeros(Nsteps-istart)
    vacf_z   = np.zeros(Nsteps-istart)
    vacf_par = np.zeros(Nsteps-istart)
    
    
    # Separated by seed:
    Nins          = []
    # This is not squared, obviously:
    alltimes  = []
    
    weirdskips   = 0 
    skippedfiles = 0
    
    for confignr in confignrs:
        print('On config number:', confignr)
        infilename_config = configfolder + 'data.config'+str(confignr)
        
        #print('infilename_all:',infilename_all)
        
        # Read in:
        #### Automatic part
        ## Find the extent of the polymers: Max z-coord of beads in the chains
        try:
            infile_config = open(infilename_config, "r")
        except:
            print('Oh, data-file! Where art thou?')
            print('file name that failed:', infilename_config)
            skippedfiles += 1
            continue # Skipping this file if it does not exist
        # Moving on, if the file exists
        lines = infile_config.readlines() # This takes some time
        # Getting the number of lines, etc.
        totlines = len(lines)         # Total number of lines
        lineend = totlines-1          # Index of last element
        
        # Extracting the number of atoms:
        words = lines[2].split()
        Nall = int(words[0])
        N    = Nall
        
        i           = linestart_data
        
        # For now: Only find the largest z-position among the beads in the chain. # Otherwise I must read off the atom number (will do that with the free bead anyways.)
        maxz = -1000 # No z-position is this small.
        
        time_start = time.process_time()
        counter = 0
        while i<totlines:
            words = lines[i].split()
            if len(words)!=0: # Some double testing going on...
                # Find properties
                # Order:  id  type mol ux  uy  uz  vx  vy   vz
                #         [0] [1]  [2] [3] [4] [5] [6] [7]  [8]
                ind      = int(words[0])-1 # Atom ids go from zero to N-1.
                atomtype = int(words[2]) 
                #molID    = int(words[2])
                z        = float(words[6])
                if atomtype==2: # Moving polymer bead. Test if this is larger than maxz:
                    if z>maxz:
                        maxz = z
                i+=1
            else:
                break
        extent_polymers = maxz
        maxz_av += maxz
        infile_config.close()
        
        for beadplacement in beadplacements:
            ## Find the position of the free bead: # I reuse quite a bit of code here...
            if long==True:
                infilename_free = endlocation_in+'long/'+'freeatom_confignr'+str(confignr)+'_beadplacement'+str(beadplacement)+'_long.lammpstrj'
            else:
                infilename_free = endlocation_in+'freeatom_confignr'+str(confignr)+'_beadplacement'+str(beadplacement)+'.lammpstrj'
            
            try:
                infile_free = open(infilename_free, "r")
            except:
                print('Oh, data-file! Where art thou?')
                print('file name that failed:', infilename_config)
                skippedfiles += 1
                continue # Skipping this file if it does not exist
            filescounter += 1
            
            # Moving on, if the file exists
            lines = infile_free.readlines() # This takes some time
            # Getting the number of lines, etc.
            totlines = len(lines)         # Total number of lines
            lineend = totlines-1          # Index of last element
            if totlines<10*Nsteps: # Hacky solution to a problem that should not even be here...
                weirdskips   += 1 # Need to track how often this happens
                skippedfiles += 1
                break        
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
            times     = []
            
            maxz = -1000 # No z-position is this small.
            
            time_start = time.process_time()
            counter = 0
            zs_fortesting = []
            while i<totlines:
                words = lines[i].split()
                if (words[0]=='ITEM:'):
                    if words[1]=='TIMESTEP':
                        words2 = lines[i+1].split() # The time step is on the next line
                        t = float(words2[0])
                        times.append(t)
                        i+=skiplines
                    elif words[1]=='NUMBER': # These will never kick in. 
                        i+=1
                    elif words[1]=='BOX':
                        i+=1
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
            times_single = np.arange(Nsteps)#*dt
            times_single_real = np.arange(Nsteps)*dt
            
            time_end = time.process_time()
            
            if max(zs_fortesting)>zhigh:
                continue
            if min(zs_fortesting)<zlow:
                continue
            
            Nin   = Nsteps
            alltimes.append(0)
            vx0 = vx[istart]
            vy0 = vy[istart]
            vz0 = vz[istart]
            vtot0 = np.array([vx0,vy0,vz0])
            for j in range(istart,Nin):
                i = j-istart      
                # Velocity:
                vxi   = vx[i]
                vyi   = vy[i]
                vzi   = vz[i]
                vtoti = np.array([vxi,vyi,vzi])
                vacf_tot[i] += np.dot(vtoti,vtot0)
                vacf_x[i]   += vxi*vx0
                vacf_y[i]   += vyi*vy0
                vacf_z[i]   += vzi*vz0
                vacf_par[i] += vxi*vx0+vyi*vy0
                average_counter[i] +=1
            
            Nins.append(Nin)
        
    print('filescounter:', filescounter)
    maxz_av /= filescounter
    
    print('maxz_av:', maxz_av)
    
    allRs      = np.array(allRs)
    alltimes   = np.array(alltimes)
    Nins       = np.array(Nins)
    
    Ninbrush = Nsteps # Default, in case it does not exit the brush
    for i in range(Nsteps-istart):
        counter = average_counter[i]
        if counter!=0:
            vacf_tot[i] /=counter
            vacf_x[i]   /=counter
            vacf_y[i]   /=counter
            vacf_z[i]   /=counter
            vacf_par[i] /=counter
        else:
           Ninbrush = i-1
           break
    
    times_single      = times_single[istart:Ninbrush]
    times_single_real = times_single_real[istart:Ninbrush]
    Ninbrush = Ninbrush-istart
    # Velocity autocorrelation function:
    vacf_tot = vacf_tot[0:Ninbrush]
    vacf_x   = vacf_x[0:Ninbrush]
    vacf_y   = vacf_y[0:Ninbrush]
    vacf_z   = vacf_z[0:Ninbrush]
    vacf_par = vacf_par[0:Ninbrush]
    
    ## Average, SI units:
    # Velocity:
    vacf_tot_SI  = vacf_tot*unitlength/unittime ##/timestepsize
    vacf_x_SI = vacf_x*unitlength/unittime ##/timestepsize
    vacf_y_SI = vacf_y*unitlength/unittime ##/timestepsize
    vacf_z_SI = vacf_z*unitlength/unittime ##/timestepsize
    vacf_par_SI = vacf_par*unitlength/unittime ##/timestepsize

    outfile_vacf = open(outfilename_vacf,'w')
    outfile_vacf.write('Time step; Time_SI; vacf_tot_SI;  vacf_z_SI;  vacf_par_SI;  vacf_x_SI;  vacf_y_SI\n')
    for i in range(Ninbrush):
        outfile_vacf.write('%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n' % (times_single[i],times_single_real[i],vacf_tot_SI[i],vacf_z_SI[i],vacf_par_SI[i],vacf_x_SI[i],vacf_y_SI[i])) 
    outfile_vacf.close()
    
    # Find D by the integral:
    Dtot = 0
    Dx   = 0
    Dy   = 0
    Dz   = 0
    Dpar = 0
    for i in range(Ninbrush):
        Dtot += vacf_tot_SI[i]
        Dx   += vacf_x_SI[i]
        Dy   += vacf_y_SI[i]
        Dz   += vacf_z_SI[i]
        Dpar += vacf_par_SI[i]
    Dtot *= dt
    Dx   *= dt
    Dy   *= dt
    Dz   *= dt
    Dpar *= dt
    
    print('--------------')
    print('spacing:',spacing)
    print('Dtot:',Dtot)
    print('Dx:',Dx)
    print('Dy:',Dy)
    print('Dz:',Dz)
    print('Dpar:',Dpar)

    outfile_Dvacf = open(outfilename_Dvacf,'w')
    outfile_Dvacf.write('d   Dtot   Dz   Dpar   Dx   Dy\n')
    outfile_Dvacf.write('%.2f %.5e %.5e %.5e %.5e %.5e\n' % (spacing,Dtot,Dz,Dpar,Dx,Dy))
    outfile_Dvacf.close()
    
    plt.figure(figsize=(6,5))
    plt.plot(times_single_real, vacf_tot_SI, label=r'$<\vec{v}(0)\cdot\vec{v}(t)>$')
    plt.plot(times_single_real, vacf_x_SI, label=r'$<v_x(0)\cdot v_x(t)>$')
    plt.plot(times_single_real, vacf_y_SI, label=r'$<v_y(0)\cdot v_y(t)>$')
    plt.plot(times_single_real, vacf_z_SI, label=r'$<v_z(0)\cdot v_z(t)>$')
    plt.plot(times_single_real, vacf_par_SI, label=r'$<v_\parallel\cdot v_\parallel(t)>$')
    plt.xlabel(r'Time [s]')
    plt.ylabel(r'VACF')
    plt.title(r'VACF vs time, d=%.2f, start index %i' % (spacing,starti))
    plt.tight_layout()
    plt.savefig(plotname)