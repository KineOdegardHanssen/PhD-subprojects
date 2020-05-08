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
# Settings
Nconfigs    = 100#100
Nplacements = 10#10
#
damp = 10
# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
popup_plots = False
spacing = 5
psigma  = 1#.5
density = 0.238732414637843 # Yields mass 1 for bead of radius 1 nm
#pmass   = 1.5
print('spacing:', spacing)
print('psigma:', psigma)
# I need to set the file name in an easier way, but for now I just use this:  ## Might want to add a loop too, if I have more files...

# Should divide into folders in a more thorough manner?
# Extracting the correct names (all of them)
T        = 3
plotseed = 0
plotdirs = False
test_sectioned = False
#seeds  = [23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]
confignrs      = np.arange(1,Nconfigs+1)#300)#20)#101)#1001)#22)
beadplacements = np.arange(1,Nplacements+1)
Nconfigs       = len(confignrs)      # So that I don't have to change that much
Nplacements    = len(beadplacements)
Nfiles         = Nconfigs*Nplacements
maxz_av        = 0
filescounter   = 0
Nsteps         = 2001 # 20001 before
#writeevery   = 10                  # I'm writing to file every this many time steps
unitlength     = 1e-9
unittime       = 2.38e-11 # s
timestepsize   = 0.00045*unittime#*writeevery
Npartitions    = 5 # For extracting more walks from one file (but is it really such a random walk here...?)
minlength      = int(floor(Nsteps/Npartitions)) # For sectioning the data
gatheredconfignrs = np.arange(1,Nfiles+1) # For some old printing
print('timestepsize:', timestepsize)
## Weird cutoff (bead):
#namebase    = '_quadr_M9N101_ljunits_spacing%i_Langevin_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_theta0is180_pmass1.5_sect_placeexact_ljcut1p122' %spacing
#folderbase  = 'Part_in_chgr_subst_all_quadr_M9N101_ljunits_Langevin_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_theta0is180_pmass1.5_sect_placeexact_ljcut1p122'
# Usual cutoff (bead):
#foldername  = 'Quadr_M9N101_ljunits_Langevin_scaled_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_ljcut1p122'
#endlocation = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Brush/'+foldername+'/Spacing'+str(spacing)+'/damp50_diffseedLgv'+'/Sigma_bead_'+str(psigma)+'/'
basepath        = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/Spacing%i/' % spacing
location_config = basepath + 'Initial_configs/Before_bead/'
endlocation     = basepath + 'Results/'
#namebase_start = '_quadr_M9N101_ljunits_'
#folderbase_mid = 'Langevin_scaled_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_pmass1.5'
#folderbase_end = '_ljcut1p122'
#namebase = namebase_start+'spacing%i_' % spacing + folderbase_mid + '_psigma' +str(psigma)+'_sect_placeexact_ljcut1p122'
#namebase = 'spacing%i_Langevin_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_pmass1.5_psigma' % spacing +str(psigma)+'_pstdcutoff_sect_placeexact_ljcut1p122'
#folderbase = 'Part_in_chgr_subst'+namebase_start+folderbase_mid+folderbase_end

#endlocation       = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Brush/'+foldername+'/Spacing'+str(spacing)+'/Sigma_bead_'+str(psigma)
filestext            = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])+'_placements'+str(beadplacements[0])+'to'+str(beadplacements[-1])
#outfilename          = endlocation+'lammpsdiffusion_qdrgr_'+namebase+filestext+'.txt' # This is redundant
# Text files
outfilename_ds       = endlocation+'av_ds_'+filestext+'.txt'                        #'lammpsdiffusion_qdrgr_'+namebase+filestext+'_av_ds.txt'
outfilename_gamma    = endlocation+'zimportance_'+filestext+'.txt'                  #'lammpsdiffusion_qdrgr_'+namebase+filestext+'_zimportance.txt'
outfilename_sections = endlocation+'sections_'+filestext+'.txt'                     #'lammpsdiffusion_qdrgr_'+namebase+filestext+'_sections.txt'
outfilename_maxz     = endlocation+'maxz_az_'+filestext+'.txt'
outfilename_dt       = endlocation+'dt.txt'

# Plots
plotname             = endlocation+filestext+'.png'                                 #'lammpsdiffusion_qdrgr_'+namebase+filestext+'.png'
plotname_all         = endlocation+'all_'+filestext+'.png'                          #'lammpsdiffusion_qdrgr_'+namebase+filestext+'_all.png'
plotname_gamma       = endlocation+'zimportance_'+filestext+'.png'                  #'lammpsdiffusion_qdrgr_'+namebase+filestext+'_zimportance.png'
plotname_SI          = endlocation+'SI_'+filestext+'.png'                           #'lammpsdiffusion_qdrgr_'+namebase+filestext+'_SI.png'
plotname_parallel_orthogonal = endlocation+'par_ort_'+filestext+'.png'              #'lammpsdiffusion_qdrgr_'+namebase+filestext+'_par_ort.png' 
plotname_dx_dy_dz    = endlocation+'dx_dy_dz_'+filestext+'.png'                     #'lammpsdiffusion_qdrgr_'+namebase+filestext+'_dx_dy_dz.png' 
plotname_parallel    = endlocation+'par_'+filestext+'.png'                          #'lammpsdiffusion_qdrgr_'+namebase+filestext+'_par.png' 
plotname_orthogonal  = endlocation+'ort_'+filestext+'.png'                          #'lammpsdiffusion_qdrgr_'+namebase+filestext+'_ort.png'
plotname_short_all   = endlocation+'short_all_'+filestext+'.png'                    #'lammpsdiffusion_qdrgr_'+namebase+filestext+'_short_all.png'
plotname_velocity    = endlocation+'velocity_'+filestext+'.png'                     #'lammpsdiffusion_qdrgr_'+namebase+filestext+'_velocity.png'
plotname_velocity_SI = endlocation+'velocity_SI_'+filestext+'.png'                  #'lammpsdiffusion_qdrgr_'+namebase+filestext+'_velocity_SI.png'
plotname_velocity_sq = endlocation+'velocity_sq_'+filestext+'.png'                  #'lammpsdiffusion_qdrgr_'+namebase+filestext+'_velocity_sq.png'
plotname_sectioned_average = endlocation+'sections_'+filestext+'.png'               #'lammpsdiffusion_'+namebase+filestext+'_sections.png'
plotname_sectioned_average_vs_steps = endlocation+'sections_steps_'+filestext+'.png'#'lammpsdiffusion_'+namebase+filestext+'_sections_steps.png'
plotname_traj_xy = endlocation+'traj_xy_config'+str(confignrs[-1])+'.png'           #'lammpsdiffusion_'+namebase+filestext+'_traj_xy.png'
plotname_traj_xz = endlocation+'traj_xz_config'+str(confignrs[-1])+'.png'           #'lammpsdiffusion_'+namebase+filestext+'_traj_xz.png'
plotname_traj_yz = endlocation+'traj_yz_config'+str(confignrs[-1])+'.png'           #'lammpsdiffusion_'+namebase+filestext+'_traj_yz.png'
plotname_traj_xt = endlocation+'traj_xt_config'+str(confignrs[-1])+'.png'           #'lammpsdiffusion_'+namebase+filestext+'_traj_xt.png'
plotname_traj_yt = endlocation+'traj_yt_config'+str(confignrs[-1])+'.png'           #'lammpsdiffusion_'+namebase+filestext+'_traj_yt.png'
plotname_traj_zt = endlocation+'traj_zt_config'+str(confignrs[-1])+'.png'           #'lammpsdiffusion_'+namebase+filestext+'_traj_zt.png'
plotname_th_hist = endlocation+'th_hist_'+filestext+'.png'                          #'lammpsdiffusion_'+namebase+filestext+'_th_hist.png'


## Setting arrays
# Prepare for sectioning distance data:
time_walks_SI, steps, partition_walks, numberofsamples, len_all, lengths, startpoints = datr.partition_holders_averaged(Nsteps,minlength)

print('Sections, walks:')
print('minlength:', minlength)
print('lenghts:', lengths)
# These are all squared:
# All together:
allRs     = []
alldxs    = []
alldys    = []
alldzs    = []
# Averages:
averageR2s  = np.zeros(Nsteps)   # Distances, squared
averagedx2s = np.zeros(Nsteps)
averagedy2s = np.zeros(Nsteps)
averagedz2s = np.zeros(Nsteps)
averagedxs  = np.zeros(Nsteps)   # Distances, non-squared
averagedys  = np.zeros(Nsteps)
averagedzs  = np.zeros(Nsteps)
averagev2s  = np.zeros(Nsteps)  # Velocities, squared
averagevx2s = np.zeros(Nsteps)
averagevy2s = np.zeros(Nsteps)
averagevz2s = np.zeros(Nsteps)
averagevs   = np.zeros(Nsteps)  # Velocities
averagevxs  = np.zeros(Nsteps)
averagevys  = np.zeros(Nsteps)
averagevzs  = np.zeros(Nsteps)
averagevparallel  = np.zeros(Nsteps)
averagevparallel2 = np.zeros(Nsteps) # Velocities, squared
averagedparallel2 = np.zeros(Nsteps) # Distances, squared
average_counter   = np.zeros(Nsteps)
average_walks     = copy.deepcopy(partition_walks)
average_walks_SI  = copy.deepcopy(partition_walks)
average_counters  = copy.deepcopy(partition_walks) # Okay, this name is confusing.


# Separated by seed:
Rs_byseed    = []
dxs_byseed   = []
dys_byseed   = []
dzs_byseed   = []
gamma_byseed = []
times_byseed = []
gamma_avgs   = []
single_slopes = []
rmsds         = []
Nins          = []
#
sections_walk  = []
sections_steps = []
# This is not squared, obviously:
alltimes  = []

linestart_data = 22
skippedfiles = 0

for confignr in confignrs:
    print('On config number:', confignr)
    infilename_config  = location_config+'data.config'+str(confignr)
    plotname_dirs      = endlocation+'dxdydzR2_config'+str(confignr)+'.png'
    plotname_testsect  = endlocation+'testsectioned_config'+str(confignr)+'.png'
    
    
    # Read in:
    #### Automatic part
    ## Find the extent of the polymers: Max z-coord of beads in the chains
    try:
        infile_config = open(infilename_config, "r")
    except:
        print('Oh, data-file! Where art thou?')
        print('file name that failed:', infilename_config)
        skippedfiles += Nplacements
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
        #print('beadplacement:', beadplacement, ' of', Nplacements)
        ## Find the position of the free bead:
        infilename_free = basepath+'freeatom_confignr'+str(confignr)+'_beadplacement'+str(beadplacement)+'.lammpstrj' #endlocation+'freeatom_density'+str(density)+'_confignr'+str(confignr)+'.lammpstrj'
        
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
        while i<totlines:
            words = lines[i].split()
            if (words[0]=='ITEM:' and words[1]=='TIMESTEP'): # Some double testing going on...
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
                counter+=1
                i+=1
        infile_free.close()
        dt = (times[1]-times[0])*timestepsize # This might be handy
        #print('dt:', dt)
        #print('times[1]-times[0]:',times[1]-times[0])
        times_single = np.arange(Nsteps)#*dt
        times_single_real = np.arange(Nsteps)*dt
        
        time_end = time.process_time()
        
        pos_inpolymer = []
        
        #before_in
        
        # I do not take into account that the free bead can enter the polymer grid and then exit again. If that happens, there might be discontinuities or weird kinks in the data (if you look at the graphs...) # I do not have this in my test-dataset. Maybe make the bead lighter and see what happens?
         
        Nin   = 0
        maxzpol = 0
        for i in range(counter):
            thesepos = positions[i]
            z        = thesepos[2]
            if z>extent_polymers: # If the polymer is in bulk # We don't want it to go back and forth between brush and bulk # That will cause discontinuities in our data
                break
            else:
                pos_inpolymer.append(thesepos)
                if z>maxzpol:
                    maxzpol = z
                Nin+=1
        
        startpos_in   = pos_inpolymer[0]
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
        gamma_temp = []
        
        R_temp.append(0)
        dx_temp.append(0)
        dy_temp.append(0)
        dz_temp.append(0)
        step_temp.append(0)
        allRs.append(0)        # We will set this here since we know the value
        alltimes.append(0)
        for i in range(1,Nin):       
            this_in = pos_inpolymer[i]
            dist = this_in-startpos_in
            dx   = dist[0]         # Signed
            dy   = dist[1]
            dz   = dist[2]
            R2   = np.dot(dist,dist) # Squared
            dx2  = dx*dx
            dy2  = dy*dy
            dz2  = dz*dz
            gamma = (R2-dz2)/R2
            # Velocity:
            vxi  = vx[i]
            vyi  = vy[i]
            vzi  = vz[i]
            vx2i = vxi*vxi # squared
            vy2i = vyi*vyi
            vz2i = vzi*vzi
            v2i  = vx2i + vy2i + vz2i
            # All together:
            allRs.append(R2)
            alldxs.append(dx2)
            alldys.append(dy2)
            alldzs.append(dz2)
            alltimes.append(i)
            # Averages:
            averageR2s[i] +=R2   # Distance
            averagedx2s[i]+=dx2
            averagedy2s[i]+=dy2
            averagedz2s[i]+=dz2
            averagedxs[i] +=dx   # Distance, signed
            averagedys[i] +=dy
            averagedzs[i] +=dz
            averagevs[i]  += np.sqrt(vx2i+vy2i+vz2i)  # Velocity
            averagevxs[i] += vxi
            averagevys[i] += vyi
            averagevzs[i] += vzi
            averagev2s[i] += vx2i+vy2i+vz2i            # Velocity, squared
            averagevx2s[i]+= vx2i
            averagevy2s[i]+= vy2i
            averagevz2s[i]+= vz2i
            averagevparallel[i]  += np.sqrt(vx2i+vy2i)
            averagevparallel2[i] += vx2i+vy2i
            averagedparallel2[i] += dx2+dy2 # Distance
            average_counter[i] +=1
            # Separated by seed:
            R_temp.append(R2)
            dx_temp.append(dx2)
            dy_temp.append(dy2)
            dz_temp.append(dz2)
            gamma_temp.append(gamma)
            step_temp.append(i)
        Rs_byseed.append(R_temp)   # Don't need these as of yet. Might use them later.
        dxs_byseed.append(dx_temp)
        dys_byseed.append(dy_temp)
        dzs_byseed.append(dz_temp)
        gamma_byseed.append(gamma_temp)
        times_byseed.append(step_temp)
        gamma_avg = np.mean(gamma_temp)
        gamma_avgs.append(gamma_avg)
        
        coeffs = np.polyfit(step_temp,R_temp,1)
        a = coeffs[0]
        b = coeffs[1]
        
        linefit = a*np.array(step_temp)+b
        
        rmsd_with_linefit = rmsd(R_temp,linefit)
        
        single_slopes.append(a)
        rmsds.append(rmsd_with_linefit)
        Nins.append(Nin)
        
        if (Nfiles<11 and plotdirs==True): # Don't plot too many files
            plt.figure(figsize=(6,5))
            plt.plot(step_temp, R_temp, label=r'$<R^2>$')
            plt.plot(step_temp, dx_temp, label=r'$<dx^2>$')
            plt.plot(step_temp, dy_temp, label=r'$<dy^2>$')
            plt.plot(step_temp, dz_temp, label=r'$<dz^2>$')
            plt.xlabel(r'Step number')
            plt.ylabel(r'Distance$^2$ [in unit length]')
            plt.title('Random walk in polymer brush, config %i, d = %i nm' % (confignr,spacing))
            plt.tight_layout()
            plt.legend(loc='upper left')
            plt.savefig(plotname_dirs)
        
        # Partitioned:
        # Npartitions
        part_walks    = []
        part_steps    = []
        part_start    = 0
        #smallest_part_walk = int(math.floor(Nin/Npartitions))
        for i in range(Npartitions):
            this_part    = []
            these_steps  = []
            part_start   = startpoints[i]
            len_thiswalk = lengths[i]
            part_end     = part_start+len_thiswalk
            if part_start>(Nin-1):              # Do not start this section if the bead has left the brush.
                break
            rstart = pos_inpolymer[part_start]
            for j in range(len_thiswalk):
                #print('i:',i,', j:',j)
                if part_start+j>(Nin-1):
                    break
                rthis =  pos_inpolymer[part_start+j]
                drvec = rthis - rstart
                dr2   = np.dot(drvec,drvec)
                '''
                print('part_start:',part_start, '; part_start+j:', part_start+j, ';   j:',j)
                print('rthis:',rthis)
                print('rstart:',rstart)
                print('drvec:', drvec)
                print('dr2:', dr2)
                '''
                average_walks[i][j]    +=dr2#= average_walks[i][j] + dr2 # Notation? [i,j] or [i][j] # Ugh, and should I have one of these for R2, dx2, dy2 and dz2?
                average_counters[i][j] +=1#= average_counters[i][j] + 1
                this_part.append(dr2)
                these_steps.append(j)
            part_walks.append(this_part)
            part_steps.append(these_steps)
        sections_walk.append(part_walks)
        sections_steps.append(part_steps)
        
        if test_sectioned==True and config==1: # Do not make too many plot
            plt.figure(figsize=(6,5))
            for i in range(Npartitions):
                plt.plot(part_steps[i], part_walks[i], label='Section %i' % i)
            plt.xlabel(r'Step number')
            plt.ylabel(r'Distance$^2$ [in unit length]')
            plt.title('Random walk in polymer brush, config %i, d = %i nm' % (confignr,spacing))
            plt.tight_layout()
            plt.legend(loc='upper left')
            plt.savefig(plotname_testsect)
        # Do something else than this? Use the data I already have?
        #time_beforepartition = time.process_time()
        #partition_walks = datr.partition_dist_averaged(positions, partition_walks, Nsteps, minlength) # partition_walks is being updated inside of the function
        #time_afterpartition = time.process_time()
        #print('Time, partitioning:',time_afterpartition-time_beforepartition)
    
print('filescounter:', filescounter)
maxz_av /= filescounter

newlen = len(gamma_avgs)

allRs      = np.array(allRs)
alltimes   = np.array(alltimes)
gamma_avgs = np.array(gamma_avgs)
single_slopes = np.array(single_slopes)
rmsds         = np.array(rmsds)
Nins          = np.array(Nins)

# Sorted arrays:
sorted_gamma_avgs = np.sort(gamma_avgs)
plotgammaagainst  = np.arange(newlen)
index_sorted      = np.argsort(gamma_avgs)

# Making arrays
rmsds_sorted         = np.zeros(newlen)
seeds_sorted         = np.zeros(newlen)
Nins_sorted          = np.zeros(newlen)
single_slopes_sorted = np.zeros(newlen)

# Calculating averages:
#Know the velocity at t=0, same for all
vxi = vx[0]
vyi = vy[0]
vzi = vz[0]
averagevs[0]  = np.sqrt(vxi*vxi + vyi*vyi + vzi*vzi)
averagevxs[0] = vxi
averagevys[0] = vyi
averagevzs[0] = vzi
averagevparallel[0] = np.sqrt(vxi*vxi+vyi*vyi)

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
        # Velocity
        averagevs[i]/=counter
        averagevxs[i]/=counter
        averagevys[i]/=counter
        averagevzs[i]/=counter
        averagevparallel[i]/=counter
        # Velocity squared
        averagev2s[i]/=counter
        averagevx2s[i]/=counter
        averagevy2s[i]/=counter
        averagevz2s[i]/=counter
        averagevparallel2[i]/=counter
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
# Velocity to the first power
averagevs  = averagevs[0:Ninbrush]
averagevxs = averagevxs[0:Ninbrush]
averagevys = averagevys[0:Ninbrush]
averagevzs = averagevzs[0:Ninbrush]
averagevparallel = averagevparallel[0:Ninbrush]
# Velocity squared
averagev2s  = averagev2s[0:Ninbrush]
averagevx2s = averagevx2s[0:Ninbrush]
averagevy2s = averagevy2s[0:Ninbrush]
averagevz2s = averagevz2s[0:Ninbrush]
averagevparallel2 = averagevparallel2[0:Ninbrush]

print('len(averageR2s):',len(averageR2s))
print('Nsteps:', Nsteps)

# Sectioned walks
for i in range(Npartitions):
    for j in range(lengths[i]):
        if average_counters[i][j]!=0:
            #print('Before division: Average_walks=', average_walks[i][j], '     average_counters=',average_counters[i][j])
            average_walks[i][j]/=average_counters[i][j]
            #print('After division: Average_walks=', average_walks[i][j], '     average_counters=',average_counters[i][j])

if Nfiles<12:
    print('gamma_avgs:',gamma_avgs)
    print('index_sorted:',index_sorted)
# Opening file
outfile_gamma = open(outfilename_gamma, 'w')
outfile_gamma.write('This is an attempt to find the relative contribution of dz to R. The smaller the value in the second column, the more important dz is.\nOrder: <(R^2(n)-dz^2(n))/R^2(n)>_n, slope a of line fit, rmsd R^2(n) line fit, Nin, confignr\n')
for i in range(newlen):
    index = index_sorted[i]
    rmsds_sorted[i] = rmsds[index]
    seeds_sorted[i] = gatheredconfignrs[index]
    Nins_sorted[i]  = Nins[index]
    single_slopes_sorted[i] = single_slopes[index]
    outfile_gamma.write('%.16f %.16f %.16f %i %i\n' % (sorted_gamma_avgs[i], single_slopes_sorted[i], rmsds_sorted[i], Nins_sorted[i], seeds_sorted[i]))
outfile_gamma.close()


## Finding the diffusion coefficient (in polymer)
# Finding the last 25% of the data:               # Do this later! In loop! # But do now just for testing purposes.
# Performing the line fit:
coeffs_poly, covs = polyfit(alltimes, allRs, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a_poly = coeffs_poly[0]
b_poly = coeffs_poly[1]
rms_D_poly = np.sqrt(covs[0,0])/6.
rms_b_poly = np.sqrt(covs[1,1])
D_poly = a_poly/6.

'''
coeffs, covs = np.polyfit(insteps, indata, 1, full=False, cov=True) # Using all the data in the fit
    av_slope[j]  = coeffs[0]
    rms_slope[j] = np.sqrt(covs[0,0])   
'''
 
fit_poly = a_poly*alltimes+b_poly

# Line fit to averaged data:
coeffs_poly, covs = polyfit(times_single, averageR2s, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a_poly_av = coeffs_poly[0]
b_poly_av = coeffs_poly[1]
rms_D_poly_av = np.sqrt(covs[0,0])/6.
rms_b_poly_av = np.sqrt(covs[1,1])
D_poly_av = a_poly_av/6.
 
fit_poly_av = a_poly_av*times_single+b_poly_av


print('Fit, all data:')
print('b_poly (should be small):', b_poly)
print('D_poly:', D_poly)

print('Fit, average:')
print('b_poly (should be small):', b_poly_av)
print('D_poly:', D_poly_av)

## Average, SI units:
# Distance:
averageR2s_SI  = averageR2s*unitlength**2
averagedx2s_SI = averagedx2s*unitlength**2
averagedy2s_SI = averagedy2s*unitlength**2
averagedz2s_SI = averagedz2s*unitlength**2
averagedparallel2_SI = averagedparallel2*unitlength**2
# Velocity:
averagevs_SI  = averagevs*unitlength/timestepsize
averagevxs_SI = averagevxs*unitlength/timestepsize
averagevys_SI = averagevys*unitlength/timestepsize
averagevzs_SI = averagevzs*unitlength/timestepsize
averagevparallel_SI = averagevparallel*unitlength/timestepsize

print('----------------------------------')
print(' averageR2s_SI[2]:', averageR2s_SI[2])

coeffs_poly, covs = polyfit(times_single_real, averageR2s_SI, 1, full=False, cov=True) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
a_poly_SI = coeffs_poly[0]
b_poly_SI = coeffs_poly[1]
rms_D_poly_SI = np.sqrt(covs[0,0])/6.
rms_b_poly_SI = np.sqrt(covs[1,1])
D_poly_SI = a_poly_SI/6.

print(' averageR2s_SI[2]:', averageR2s_SI[2])

print('----------------------------------')
fit_poly_SI = a_poly_SI*times_single_real+b_poly_SI

print('Fit, average, SI:')
print('b_poly (should be small):', b_poly_SI)
print('D_poly:', D_poly_SI)

# Do not perform fit here:
'''
outfile = open(outfilename, 'w')
outfile.write('D_poly: %.16f %.16f\n' % (D_poly,rms_D_poly))
outfile.write('b_poly: %.16f %.16f\n' % (b_poly,rms_b_poly))
outfile.close()
'''

outfile_ds = open(outfilename_ds,'w')
outfile_ds.write('Time step; Time; <R^2>; <dx^2>; <dy^2>; <dz^2>; <dx^2+dy^2>\n')
for i in range(len(averageR2s)):
    outfile_ds.write('%i %.16e %.16e %.16e %.16e %.16e %.16e\n' % (times_single[i], times_single_real[i], averageR2s_SI[i], averagedx2s_SI[i], averagedy2s_SI[i], averagedz2s_SI[i], averagedparallel2_SI[i]))
outfile_ds.close()

xpos = []
ypos = []
zpos = []

outfilename_xoft     = endlocation+'xoft_config%i.txt' % confignr
outfilename_yoft     = endlocation+'yoft_config%i.txt' % confignr
outfilename_zoft     = endlocation+'zoft_config%i.txt' % confignr


outfile_xoft = open(outfilename_xoft,'w')
outfile_yoft = open(outfilename_yoft,'w')
outfile_zoft = open(outfilename_zoft,'w')

print('len(positions):',len(positions))
print('len(times_single):',len(times_single))
'''# I still have some trouble when I don't have enough files (maybe when the last one is not made yet?)
Nloop = len(positions)
if len(times_single)<Nloop:
    Nloop = len(times_single)
'''
for i in range(len(times_single)):
    pos = positions[i]
    xthis = pos[0]
    ythis = pos[1]
    zthis = pos[2]
    xpos.append(xthis)
    ypos.append(ythis)
    zpos.append(zthis)
    outfile_xoft.write('%i %.16f\n' % (times_single[i],xthis))
    outfile_yoft.write('%i %.16f\n' % (times_single[i],ythis))
    outfile_zoft.write('%i %.16f\n' % (times_single[i],zthis))
outfile_xoft.close()
outfile_yoft.close()
outfile_zoft.close()



# Trajectory in plane
plt.figure(figsize=(6,5))
plt.plot(xpos, ypos, '.')
plt.xlabel(r'x')
plt.ylabel(r'y')
plt.title('Trajectory in xy-plane')
plt.tight_layout()
plt.savefig(plotname_traj_xy)

plt.figure(figsize=(6,5))
plt.plot(xpos, zpos, '.')
plt.xlabel(r'x')
plt.ylabel(r'z')
plt.title('Trajectory in xz-plane')
plt.tight_layout()
plt.savefig(plotname_traj_xz)

plt.figure(figsize=(6,5))
plt.plot(ypos, zpos, '.')
plt.xlabel(r'y')
plt.ylabel(r'z')
plt.title('Trajectory in yz-plane')
plt.tight_layout()
plt.savefig(plotname_traj_yz)

### x, y, z vs t
plt.figure(figsize=(6,5))
plt.plot(times_single, xpos, '.')
plt.xlabel(r't')
plt.ylabel(r'x')
plt.title('x(t)')
plt.tight_layout()
plt.savefig(plotname_traj_xt)

plt.figure(figsize=(6,5))
plt.plot(times_single, ypos, '.')
plt.xlabel(r't')
plt.ylabel(r'y')
plt.title('y(t)')
plt.tight_layout()
plt.savefig(plotname_traj_yt)

plt.figure(figsize=(6,5))
plt.plot(times_single, zpos, '.')
plt.xlabel(r't')
plt.ylabel(r'z')
plt.title('z(t)')
plt.tight_layout()
plt.savefig(plotname_traj_zt)


outfile_sections = open(outfilename_sections, 'w')
for i in range(Npartitions):
    outfile_sections.write('Section %i\n' % i) # When reading from the file afterwards, I can extract the section number by if words[0]=='Section' and use that to divide the data into sections
    for j in range(lengths[i]):
        average_walks_SI[i][j] = average_walks[i][j]*unitlength 
        time_walks_SI[i][j]    = steps[i][j]*timestepsize
        #print('time_walks_SI[i][j]:', time_walks_SI[i][j], ': steps[i][j]:',steps[i][j], ';   timestepsize:',timestepsize, '; steps[i][j]*timestepsize:', steps[i][j]*timestepsize)
        outfile_sections.write('%.16e %.16e\n' % (time_walks_SI[i][j],average_walks_SI[i][j]))
outfile_sections.close()


outfile_maxz = open(outfilename_maxz, 'w')
outfile_maxz.write('%.16f' % maxz_av)         # Should I have some rms value here too?
outfile_maxz.close()

outfile_dt = open(outfilename_dt, 'w')
outfile_dt.write('%.16e' % dt)
outfile_dt.close()

if popup_plots==True:
    maxind = 100#00
    plt.figure(figsize=(6,5))
    plt.plot(times_single_real, averagevs_SI, label=r'$<R^2>$')
    plt.plot(times_single_real, averagevxs_SI, label=r'$<dx^2>$')
    plt.plot(times_single_real, averagevys_SI, label=r'$<dy^2>$')
    plt.plot(times_single_real, averagevzs_SI, label=r'$<dz^2>$')
    plt.plot(times_single_real, averagevparallel_SI, label=r'$<dx^2+dy^2>$')
    plt.xlabel(r'Index (s)')
    plt.ylabel(r'Distance$^2$ [in unit length]')
    plt.title(r'RMSD in brush, d = %i nm, $\sigma_b=%.2f$, SI' % (spacing,psigma))
    plt.tight_layout()
    plt.legend(loc='upper left')
    plt.axis([0, times_single_real[maxind], min(averagevxs_SI[0:maxind]), max(averagevs_SI[0:maxind])])
    plt.show()

'''
plt.figure()
for i in range(Npartitions):
    plt.plot(steps[:][i],average_walks[:][i], label='Start at step %i' % startpoints[i])
plt.xlabel(r'Time step')
plt.ylabel(r'Distance$^2$ [m]')
plt.title('RMSD in brush, d = %i nm, SI, sectioning' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.show()


maxind = 10000
plt.figure(figsize=(6,5))
plt.plot(times_single_real, averageR2s_SI, label=r'$<R^2>$')
plt.plot(times_single_real, averagedz2s_SI, label=r'$<dz^2>$')
plt.plot(times_single_real, averagedparallel2_SI, label=r'$<dx^2+dy^2>$')
plt.xlabel(r'Index (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %i nm, $\sigma_b=%.2f$, SI' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.axis([0, times_single_real[maxind], 0, max(averageR2s_SI[0:maxind])])
plt.show()
'''
# To determine the range
xmax_plot = 100
plt.figure(figsize=(6,5))
plt.plot(times_single, averageR2s, label=r'$<R^2>$')
plt.plot(times_single, averagedx2s, label=r'$<dx^2>$')
plt.plot(times_single, averagedy2s, label=r'$<dy^2>$')
plt.plot(times_single, averagedz2s, label=r'$<dz^2>$')
plt.plot(times_single, averagedparallel2, label=r'$<dx^2+dy^2>$')
plt.xlabel(r'Index (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.axis([0, xmax_plot, 0, max(averageR2s[0:xmax_plot])])
plt.savefig(plotname_short_all)

plt.figure(figsize=(6,5))
plt.plot(alltimes, allRs, ',', label='Data, brush')
plt.plot(times_single, fit_poly_av, '--', label='Fit, brush')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.plot(alltimes, allRs, ',', label='Data, brush')
plt.plot(times_single, averageR2s, label='Average, brush')
plt.plot(alltimes, fit_poly, '--', label='Fit, data, brush')
plt.plot(times_single, fit_poly_av, '--', label='Fit, average, brush')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_all)

plt.figure(figsize=(6,5))
plt.plot(times_single, averageR2s, label=r'$<R^2>$')
plt.plot(times_single, averagedx2s, label=r'$<dx^2>$')
plt.plot(times_single, averagedy2s, label=r'$<dy^2>$')
plt.plot(times_single, averagedz2s, label=r'$<dz^2>$')
plt.plot(times_single, averagedparallel2, label=r'$<dx^2+dy^2>$')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'Averaged RMSD in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_parallel_orthogonal)

plt.figure(figsize=(6,5))
plt.plot(times_single, averagedzs, ',')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'$<dz^2>$ in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.savefig(plotname_orthogonal)

plt.figure(figsize=(6,5))
plt.plot(times_single, averagedxs, label=r'$<dx>$')
plt.plot(times_single, averagedys, label=r'$<dy>$')
plt.plot(times_single, averagedzs, label=r'$<dz>$')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance [in unit length]')
plt.title(r'Averaged RMSD in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_dx_dy_dz)

plt.figure(figsize=(6,5))
plt.plot(times_single, averagedparallel2, ',')
plt.xlabel(r'Step number')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'$<dx^2+dy^2>$ in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.savefig(plotname_parallel)

plt.figure(figsize=(6,5))
plt.plot(plotgammaagainst, sorted_gamma_avgs)
plt.xlabel(r'Different runs')
plt.ylabel(r'$\gamma = <\frac{R^2-dz^2}{R^2}>$')
plt.title(r'Random walk in polymer brush, $\gamma$; d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.savefig(plotname_gamma)

plt.figure(figsize=(6,5))
plt.plot(times_single_real, averageR2s_SI, ',', label='Average, brush')
plt.plot(times_single_real, fit_poly_SI, '--', label='Fit, average, brush')
plt.xlabel(r'Time (s)')
plt.ylabel(r'Distance$^2$ [in unit length]')
plt.title(r'RMSD in brush, d = %i nm, $\sigma_b=%.2f$' % (spacing,psigma))
plt.tight_layout()
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.legend(loc='upper left')
plt.savefig(plotname_SI)

plt.figure()
plt.plot(times_single, averagevs, label=r'$<v>$')
plt.plot(times_single, averagevxs, label=r'$<v_x>$')
plt.plot(times_single, averagevys, label=r'$<v_y>$')
plt.plot(times_single, averagevzs, label=r'$<v_z>$')
plt.plot(times_single, averagevparallel, label=r'$<v_\parallel>$')
plt.xlabel(r'Step number')
plt.ylabel(r'Velocity [in unit length/time]')
plt.title('Averaged velocity in brush, d = %i nm' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_velocity)

plt.figure()
plt.plot(times_single, averagev2s, label=r'$<v^2>$')
plt.plot(times_single, averagevx2s, label=r'$<v_x^2>$')
plt.plot(times_single, averagevy2s, label=r'$<v_y^2>$')
plt.plot(times_single, averagevz2s, label=r'$<v_z^2>$')
plt.plot(times_single, averagevparallel2, label=r'$<v_x^2+vy^2>$')
plt.xlabel(r'Step number')
plt.ylabel(r'Velocity [(length/time)$^2$]')
plt.title('Averaged velocity in brush, d = %i nm' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig(plotname_velocity_sq)

plt.figure()
plt.plot(times_single_real, averagevs_SI, label=r'v')
plt.plot(times_single_real, averagevxs_SI, label=r'vx')
plt.plot(times_single_real, averagevys_SI, label=r'vy')
plt.plot(times_single_real, averagevzs_SI, label=r'vz')
plt.xlabel(r'Time [s]')
plt.ylabel(r'Velocity [m/s]')
plt.title('Averaged velocity in brush, d = %i nm, SI' % spacing)
plt.tight_layout()
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.legend(loc='upper left')
plt.savefig(plotname_velocity_SI)




'''
print('Original:')
print('steps[0][:]:',steps[0][:])
print('steps[:][0]:',steps[:][0])
print('steps[1][:]:',steps[1][:])
print('steps[:][1]:',steps[:][1])
print('steps[2][:]:',steps[2][:])
print('steps[:][2]:',steps[:][2])


print('SI!!!!:')
print('time_walks_SI[0][:]:',time_walks_SI[0][:])
print('time_walks_SI[:][0]:',time_walks_SI[:][0])
print('time_walks_SI[1][:]:',time_walks_SI[1][:])
print('time_walks_SI[:][1]:',time_walks_SI[:][1])
print('time_walks_SI[2][:]:',time_walks_SI[2][:])
print('time_walks_SI[:][2]:',time_walks_SI[:][2])


print('average_walks_SI[0][:]:',average_walks_SI[0][:])
print('average_walks_SI[:][0]:',average_walks_SI[:][0])
print('average_walks_SI[1][:]:',average_walks_SI[1][:])
print('average_walks_SI[:][1]:',average_walks_SI[:][1])
print('average_walks_SI[2][:]:',average_walks_SI[2][:])
print('average_walks_SI[:][2]:',average_walks_SI[:][2])
'''
print('average_walks:',average_walks)
print('average_walks_SI:',average_walks_SI)
#'''

plt.figure()
for i in range(Npartitions):
    plt.plot(time_walks_SI[:][i],average_walks_SI[:][i], label='Start at step %i' % startpoints[i])
plt.xlabel(r'Time [s]')
plt.ylabel(r'Distance$^2$ [m]')
plt.title('RMSD in brush, d = %i nm, SI, sectioning' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.savefig(plotname_sectioned_average)

plt.figure()
for i in range(Npartitions):
    plt.plot(steps[:][i],average_walks_SI[:][i], label='Start at step %i' % startpoints[i])
plt.xlabel(r'Time step')
plt.ylabel(r'Distance$^2$ [m]')
plt.title('RMSD in brush, d = %i nm, SI, sectioning' % spacing)
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.savefig(plotname_sectioned_average_vs_steps)

plotname_sectioned_average_vs_steps

print('counters:',average_counter)

print('spacing:', spacing)
print('psigma:', psigma)
