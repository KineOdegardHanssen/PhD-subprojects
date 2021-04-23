import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob

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

def find_start_positive(x):
    N = len(x)
    for i in range(N):
        if x[i]>=0:
            return i # Does this work?


# Models:

def ordered_packings(d,sigma):
    phi = 1 - np.pi*(sigma/d)
    tau = (3-phi)/2.
    return 1./tau

def hyp_rev(d,sigma):
    phi = 1 - np.pi*(sigma/d)
    tau = 2-phi
    return 1./tau

def notmonosph(d,sigma):
    phi = 1 - np.pi*(sigma/d)
    tau = 1./np.sqrt(phi) # phi**(-1/2.)
    return 1./tau

def pshimd_spherepacking(d,sigma):
    phi = 1 - np.pi*(sigma/d)
    tau = phi**(-1/3.)
    return 1./tau

def overlapping_spheres(d,sigma):
    phi = 1 - np.pi*(sigma/d)
    tau = 1-np.log(phi)/2.
    return 1./tau

def overlapping_cylinders(d,sigma):
    phi = 1 - np.pi*(sigma/d)
    tau = 1-np.log(phi)
    return 1./tau

def het_cat(d,sigma):
    phi = 1 - np.pi*(sigma/d)
    tau = phi/(1-(1-phi)**(1/3.))
    return 1./tau

def catex_resin_membr(d,sigma):
    phi = 1 - np.pi*(sigma/d)
    tau = ((2.-phi)/phi)**2
    return 1./tau

showplots = True#False

# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
rsphere  = 1.0 #0.8
#sigma    = 5.0 # Change this for the modelling part

# Fixed parameters
psigma   = 1 # For instance 
damp     = 10
# Input booleans for file selection:
bulkdiffusion = False
substrate     = False
moresigmas    = False
big           = False
bulk_cut      = False
confignrs     = np.arange(1,1001)
filestext     = '_seed'+str(confignrs[0])+'to'+str(confignrs[-1])

endlocation = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/Nocut/'
if moresigmas==True:
    endlocation = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Moresigmas/Nocut/'

basepath_base      = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/'
endlocation_static = basepath_base+'D_vs_d/Nocut/'

bulklocation = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/'

bulkfilename  = bulklocation + 'diffusion_bulk'+filestext
if bulk_cut==True:
    bulkfilename = bulkfilename +'_cut.txt'
else:
    bulkfilename = bulkfilename +'_uncut.txt'

## Files to read
brushfilename_dyn  = endlocation + 'D_vs_d_better_rms_Nintervals10.txt'
brushfilename_stat = endlocation_static + 'D_vs_d_static_better_rms_Nintervals10.txt'

## Files to write to
rmsdfilename_RW = endlocation_static+'D_rmsd_close_packing_vs_RWgeom_norefl_rsphere'+str(rsphere)+'_maxlen1.0_modelr_findsigma.txt'
rmsdfilename_forest = endlocation_static+'D_rmsd_close_packing_vs_forest_modelr_findsigma.txt'
rmsdfilename_dyn = endlocation_static+'D_rmsd_close_packing_vs_dyn_modelr_findsigma.txt'
rmsdfilename_stat = endlocation_static+'D_rmsd_close_packing_vs_stat_modelr_findsigma.txt'
if big==False:
    ### d
    plotname_d_dyn     = endlocation_static+'D_vs_d_dyn_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma.png'
    plotname_d_stat    = endlocation_static+'D_vs_d_stat_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma.png'
    plotname_d_f       = endlocation_static+'D_vs_d_f_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma.png'
    plotname_cut_d = endlocation_static+'D_vs_d_dyn_vs_stat_vs_RWgeom_cut_norefl_rsphere'+str(rsphere)+'_modelr_findsigma.png' ## Use this?
else:
    ### d 
    plotname_d_dyn     = endlocation_static+'D_vs_d_dyn_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma.png'
    plotname_d_stat    = endlocation_static+'D_vs_d_stat_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma.png'
    plotname_d_f       = endlocation_static+'D_vs_d_f_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma.png'
    plotname_cut_d = endlocation_static+'D_vs_d_dyn_vs_stat_cut_vs_RWgeom_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma.png' ## Use this?

# Dynamic sims:
brushfile_dyn = open(brushfilename_dyn, 'r')

lines = brushfile_dyn.readlines()
N_dyn = len(lines)-1

# ds
dens_dyn = np.zeros(N_dyn)
spacings_dyn = np.zeros(N_dyn)
# Ds
DRs_dyn = np.zeros(N_dyn)
Dxs_dyn = np.zeros(N_dyn)
Dys_dyn = np.zeros(N_dyn)
Dzs_dyn = np.zeros(N_dyn)
Dparallel_dyn = np.zeros(N_dyn)
# Ds, stdv
DRs_stdv_dyn = np.zeros(N_dyn)
Dxs_stdv_dyn = np.zeros(N_dyn)
Dys_stdv_dyn = np.zeros(N_dyn)
Dzs_stdv_dyn = np.zeros(N_dyn)
Dparallel_stdv_dyn = np.zeros(N_dyn)

for i in range(1,N_dyn+1):
    words = lines[i].split()
    j = i-1
    
    d = float(words[0])
    spacings_dyn[j] = d
    dens_dyn[j] = np.pi/float(d**2)
    DRs_dyn[j]  = float(words[1])
    Dzs_dyn[j]  = float(words[3])
    Dparallel_dyn[j] = float(words[5])
    # Ds, stdv
    DRs_stdv_dyn[j] = float(words[2])
    Dzs_stdv_dyn[j] = float(words[4])
    Dparallel_stdv_dyn[j] = float(words[6])
    
brushfile_dyn.close()

##########

# Static sims:
brushfile_stat = open(brushfilename_stat, 'r')

lines = brushfile_stat.readlines()
N_stat = len(lines)-1

# ds
dens_stat     = np.zeros(N_stat)
spacings_stat = np.zeros(N_stat)
# Ds
DRs_stat = np.zeros(N_stat)
Dxs_stat = np.zeros(N_stat)
Dys_stat = np.zeros(N_stat)
Dzs_stat = np.zeros(N_stat)
Dparallel_stat = np.zeros(N_stat)
# Ds, stdv
DRs_stdv_stat = np.zeros(N_stat)
Dxs_stdv_stat = np.zeros(N_stat)
Dys_stdv_stat = np.zeros(N_stat)
Dzs_stdv_stat = np.zeros(N_stat)
Dparallel_stdv_stat = np.zeros(N_stat)

for i in range(1,N_stat+1):
    words = lines[i].split()
    j = i-1
    
    d = float(words[0])
    spacings_stat[j] = d
    dens_stat[j] = np.pi/float(d**2)
    DRs_stat[j]  = float(words[1])
    Dzs_stat[j]  = float(words[3])
    Dparallel_stat[j] = float(words[5])
    # Ds, stdv
    DRs_stdv_stat[j] = float(words[2])
    Dzs_stdv_stat[j] = float(words[4])
    Dparallel_stdv_stat[j] = float(words[6])
    
brushfile_stat.close()

###
#Bulk:

bulkfile  = open(bulkfilename, 'r')
# D_R2  sigmaD_R2 b_R2 sigmab_R2; D_z2  sigmaD_z2  b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2
bulklines = bulkfile.readlines()
bulkline  = bulklines[1]
words     = bulkline.split()

# Ds
DRs_bulk = float(words[0])
Dzs_bulk = float(words[4])
Dparallel_bulk = float(words[8])

# Ds, stdv
DRs_stdv_bulk = float(words[1])
Dzs_stdv_bulk = float(words[5])
Dparallel_stdv_bulk = float(words[9])

bulkfile.close()

# Divide by bulk:
for i in range(N_stat):
    DRnew   = DRs_stat[i]/DRs_bulk
    Dznew   = Dzs_stat[i]/DRs_bulk
    Dparnew = Dparallel_stat[i]/DRs_bulk 
    # Ds, stdv
    DRs_stdv_stat[i] = abs(DRnew)*np.sqrt((DRs_stdv_stat[i]/DRs_stat[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dzs_stdv_stat[i] = abs(Dznew)*np.sqrt((Dzs_stdv_stat[i]/Dzs_stat[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dparallel_stdv_stat[i] = abs(Dparnew)*np.sqrt((Dparallel_stdv_stat[i]/Dparallel_stat[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    # Ds
    DRs_stat[i] = DRnew
    Dzs_stat[i] = Dznew
    Dparallel_stat[i] = Dparnew

for i in range(N_dyn):
    DRnew   = DRs_dyn[i]/DRs_bulk
    Dznew   = Dzs_dyn[i]/DRs_bulk
    Dparnew = Dparallel_dyn[i]/DRs_bulk 
    # Ds, stdv
    DRs_stdv_dyn[i] = abs(DRnew)*np.sqrt((DRs_stdv_dyn[i]/DRs_dyn[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dzs_stdv_dyn[i] = abs(Dznew)*np.sqrt((Dzs_stdv_dyn[i]/Dzs_dyn[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dparallel_stdv_dyn[i] = abs(Dparnew)*np.sqrt((Dparallel_stdv_dyn[i]/Dparallel_dyn[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    # Ds
    DRs_dyn[i] = DRs_dyn[i]/DRs_bulk
    Dzs_dyn[i] = Dzs_dyn[i]/DRs_bulk
    Dparallel_dyn[i] =  Dparallel_dyn[i]/DRs_bulk

## Forest (a new addition):
endlocation_forest = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_forest/D_vs_d/Nocut/'
infilename_f  = endlocation_forest+'Dd_div_Dbulk_vs_d_2_forest_uncut.txt'
infile_f = open(infilename_f,'r')
lines = infile_f.readlines()
N_f   = len(lines)-1

# ds
dens_f     = np.zeros(N_f)
spacings_f = np.zeros(N_f)
# Ds
DRs_f = np.zeros(N_f)
Dxs_f = np.zeros(N_f)
Dys_f = np.zeros(N_f)
Dzs_f = np.zeros(N_f)
Dparallel_f = np.zeros(N_f)
# Ds, stdv
DRs_stdv_f = np.zeros(N_f)
Dxs_stdv_f = np.zeros(N_f)
Dys_stdv_f = np.zeros(N_f)
Dzs_stdv_f = np.zeros(N_f)
Dparallel_stdv_f = np.zeros(N_f)

for i in range(1,N_f+1):
    words = lines[i].split()
    j = i-1
    
    d = float(words[0])
    spacings_f[j] = d
    dens_f[j] = np.pi/float(d**2)
    DRs_f[j]  = float(words[1])
    Dzs_f[j]  = float(words[3])
    Dparallel_f[j] = float(words[5])
    # Ds, stdv
    DRs_stdv_f[j] = float(words[2])
    Dzs_stdv_f[j] = float(words[4])
    Dparallel_stdv_f[j] = float(words[6])
infile_f.close()

## RW part:
Nsteps   = 1000   # Increase #Maybe for other values of hitprob
Nreal    = 10000
folder   = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Randomwalk_fixedgeometry/'

infilename_Dbulk = folder+'2D_bulk_Nsteps%i_Nreal%i_D.txt' %(Nsteps,Nreal)
infile_Ddbulk = open(infilename_Dbulk,'r')
line = infile_Ddbulk.readline()
Dbulk = float(line.split()[0])
infile_Ddbulk.close()

folder   = folder+'No_reflection/'
infilename_Ddbulk = folder+'D_vs_d_Nsteps%i_Nreal%i_rsphere'%(Nsteps,Nreal)+str(rsphere)+'_norefl.txt' 
infile_Ddbulk = open(infilename_Ddbulk,'r')
lines = infile_Ddbulk.readlines()
#print('lines:',lines)

print('infilename_Ddbulk:',infilename_Ddbulk)
phit = []
DRW  = []
dg   = []
for line in lines:
    words = line.split()
    #print('words:',words)
    #if len(words)!=0:
    #    phit.append(float(words[0]))
    #    DRW.append(float(words[1]))
    d = float(words[0])
    dg.append(d)
    phit.append(np.pi*(rsphere**2)/d**2)
    DRW.append(float(words[1])/Dbulk)
infile_Ddbulk.close()
print('phit:',phit)


''' # All take phi (the porosity) as input
ordered_packings
hyp_rev
notmonosph
pshimd_spherepacking
overlapping_spheres
overlapping_cylinders
het_cat
catex_resin_membr
'''

# Porosity
dg = np.array(dg)
spacings_rw   = dg
spacings_stat = np.array(spacings_stat)
spacings_dyn  = np.array(spacings_dyn)
spacings_f    = np.array(spacings_f)
#porosity_rw   = 1-(np.pi*(sigma/dg)**2)
#porosity_stat = 1-(np.pi*(sigma/spacings_stat)**2)
#porosity_dyn  = 1-(np.pi*(sigma/spacings_dyn)**2)
#porosity_f    = 1-(np.pi*(sigma/spacings_f)**2)

# Find start indices for positive values # Some functions cannot handle negative porosities
#startind_f    = find_start_positive(porosity_f)
#startind_rw   = find_start_positive(porosity_rw)
#startind_dyn  = find_start_positive(porosity_dyn)
#startind_stat = find_start_positive(porosity_stat)

# Making arrays based on the indices above 
# Spacing
#d_f_positive    = spacings_f[startind_f:]
#d_rw_positive   = dg[startind_rw:]
#d_dyn_positive  = spacings_dyn[startind_dyn:]     # Need THIS
#d_stat_positive = spacings_stat[startind_stat:]
# Porosity
#porosity_f_positive    = porosity_f[startind_f:]
#porosity_rw_positive   = porosity_rw[startind_rw:]
#porosity_dyn_positive  = porosity_dyn[startind_dyn:]
#porosity_stat_positive = porosity_stat[startind_stat:]
# Diffusion (to find RMSD)
#D_rw_positive   = DRW[startind_rw:]
#D_stat_positive = Dparallel_stat[startind_stat:]
#D_dyn_positive  = Dparallel_dyn[startind_dyn:]
#D_f_positive    = Dparallel_f[startind_f:]


# Ordered packing
popt, pcov           = curve_fit(ordered_packings, spacings_rw, DRW)#, maxfev=thismaxfev)
sigma_op_rw          = popt[0]
popt, pcov           = curve_fit(ordered_packings, spacings_stat, Dparallel_stat)#, maxfev=thismaxfev)
sigma_op_stat        = popt[0]
popt, pcov           = curve_fit(ordered_packings, spacings_dyn, Dparallel_dyn,p0=7)#, maxfev=thismaxfev)
sigma_op_dyn         = popt[0]
popt, pcov           = curve_fit(ordered_packings, spacings_f, Dparallel_f)#, maxfev=thismaxfev)
sigma_op_f           = popt[0]

#sigma_op_dyn = 7
# Ds
Ds_op_rw   = ordered_packings(spacings_rw,sigma_op_rw)
Ds_op_stat = ordered_packings(spacings_stat,sigma_op_stat)
Ds_op_dyn  = ordered_packings(spacings_dyn,sigma_op_dyn)
Ds_op_f    = ordered_packings(spacings_f,sigma_op_f)
# Hyperbola of revolution
popt, pcov           = curve_fit(hyp_rev, spacings_rw, DRW)#, maxfev=thismaxfev)
sigma_hr_rw          = popt[0]
popt, pcov           = curve_fit(hyp_rev, spacings_stat, Dparallel_stat)#, maxfev=thismaxfev)
sigma_hr_stat        = popt[0]
popt, pcov           = curve_fit(hyp_rev, spacings_dyn, Dparallel_dyn, p0=5)#, maxfev=thismaxfev)
sigma_hr_dyn         = popt[0]
popt, pcov           = curve_fit(hyp_rev, spacings_f, Dparallel_f)#, maxfev=thismaxfev)
sigma_hr_f           = popt[0]

#sigma_hr_dyn = 5
# Ds
Ds_hr_rw   = hyp_rev(spacings_rw,sigma_hr_rw)
Ds_hr_stat = hyp_rev(spacings_stat,sigma_hr_stat)
Ds_hr_dyn  = hyp_rev(spacings_dyn,sigma_hr_dyn)
Ds_hr_f    = hyp_rev(spacings_f,sigma_hr_f)
# Not monosized spheres
popt, pcov           = curve_fit(notmonosph, spacings_rw, DRW)#, maxfev=thismaxfev)
sigma_nm_rw          = popt[0]
popt, pcov           = curve_fit(notmonosph, spacings_stat, Dparallel_stat)#, maxfev=thismaxfev)
sigma_nm_stat        = popt[0]
popt, pcov           = curve_fit(notmonosph, spacings_dyn, Dparallel_dyn)#, maxfev=thismaxfev)
sigma_nm_dyn         = popt[0]
popt, pcov           = curve_fit(notmonosph, spacings_f, Dparallel_f)#, maxfev=thismaxfev)
sigma_nm_f           = popt[0]
# Ds
Ds_nm_rw   = notmonosph(spacings_rw,sigma_nm_rw)
Ds_nm_stat = notmonosph(spacings_stat,sigma_nm_stat)
Ds_nm_dyn  = notmonosph(spacings_dyn,sigma_nm_dyn)
Ds_nm_f    = notmonosph(spacings_f,sigma_nm_f)
# pshimd_spherepacking
popt, pcov           = curve_fit(pshimd_spherepacking, spacings_rw, DRW)#, maxfev=thismaxfev)
sigma_ps_rw          = popt[0]
popt, pcov           = curve_fit(pshimd_spherepacking, spacings_stat, Dparallel_stat)#, maxfev=thismaxfev)
sigma_ps_stat        = popt[0]
popt, pcov           = curve_fit(pshimd_spherepacking, spacings_dyn, Dparallel_dyn)#, maxfev=thismaxfev)
sigma_ps_dyn         = popt[0]
popt, pcov           = curve_fit(pshimd_spherepacking, spacings_f, Dparallel_f)#, maxfev=thismaxfev)
sigma_ps_f           = popt[0]
# Ds
Ds_ps_rw   = pshimd_spherepacking(spacings_rw,sigma_ps_rw)
Ds_ps_stat = pshimd_spherepacking(spacings_stat,sigma_ps_stat)
Ds_ps_dyn  = pshimd_spherepacking(spacings_dyn,sigma_ps_dyn)
Ds_ps_f    = pshimd_spherepacking(spacings_f,sigma_ps_f)
# Overlapping spheres
popt, pcov           = curve_fit(overlapping_spheres, spacings_rw, DRW)#, maxfev=thismaxfev)
sigma_os_rw          = popt[0]
popt, pcov           = curve_fit(overlapping_spheres, spacings_stat, Dparallel_stat)#, maxfev=thismaxfev)
sigma_os_stat        = popt[0]
popt, pcov           = curve_fit(overlapping_spheres, spacings_dyn, Dparallel_dyn)#, maxfev=thismaxfev)
sigma_os_dyn         = popt[0]
popt, pcov           = curve_fit(overlapping_spheres, spacings_f, Dparallel_f)#, maxfev=thismaxfev)
sigma_os_f           = popt[0]
# Ds
Ds_os_rw   = overlapping_spheres(spacings_rw,sigma_os_rw)
Ds_os_stat = overlapping_spheres(spacings_stat,sigma_os_stat)
Ds_os_dyn  = overlapping_spheres(spacings_dyn,sigma_os_dyn)
Ds_os_f    = overlapping_spheres(spacings_f,sigma_os_f)
# Overlapping cylinders
popt, pcov           = curve_fit(overlapping_cylinders, spacings_rw, DRW)#, maxfev=thismaxfev)
sigma_oc_rw          = popt[0]
popt, pcov           = curve_fit(overlapping_cylinders, spacings_stat, Dparallel_stat)#, maxfev=thismaxfev)
sigma_oc_stat        = popt[0]
popt, pcov           = curve_fit(overlapping_cylinders, spacings_dyn, Dparallel_dyn)#, maxfev=thismaxfev)
sigma_oc_dyn         = popt[0]
popt, pcov           = curve_fit(overlapping_cylinders, spacings_f, Dparallel_f)#, maxfev=thismaxfev)
sigma_oc_f           = popt[0]
# Ds
Ds_oc_rw   = overlapping_cylinders(spacings_rw,sigma_oc_rw)
Ds_oc_stat = overlapping_cylinders(spacings_stat,sigma_oc_stat)
Ds_oc_dyn  = overlapping_cylinders(spacings_dyn,sigma_oc_dyn)
Ds_oc_f    = overlapping_cylinders(spacings_f,sigma_oc_f)
# Heterogeneous catalyst
popt, pcov           = curve_fit(het_cat, spacings_rw, DRW)#, maxfev=thismaxfev)
sigma_hc_rw          = popt[0]
popt, pcov           = curve_fit(het_cat, spacings_stat, Dparallel_stat)#, maxfev=thismaxfev)
sigma_hc_stat        = popt[0]
popt, pcov           = curve_fit(het_cat, spacings_dyn, Dparallel_dyn)#, maxfev=thismaxfev)
sigma_hc_dyn         = popt[0]
popt, pcov           = curve_fit(het_cat, spacings_f, Dparallel_f)#, maxfev=thismaxfev)
sigma_hc_f           = popt[0]
# Ds
Ds_hc_rw   = het_cat(spacings_rw,sigma_hc_rw)
Ds_hc_stat = het_cat(spacings_stat,sigma_hc_stat)
Ds_hc_dyn  = het_cat(spacings_dyn,sigma_hc_dyn)
Ds_hc_f    = het_cat(spacings_f,sigma_hc_f)
# Cation exchange resin memrane
popt, pcov           = curve_fit(catex_resin_membr, spacings_rw, DRW)#, maxfev=thismaxfev)
sigma_rm_rw          = popt[0]
popt, pcov           = curve_fit(catex_resin_membr, spacings_stat, Dparallel_stat)#, maxfev=thismaxfev)
sigma_rm_stat        = popt[0]
popt, pcov           = curve_fit(catex_resin_membr, spacings_dyn, Dparallel_dyn)#, maxfev=thismaxfev)
sigma_rm_dyn         = popt[0]
popt, pcov           = curve_fit(catex_resin_membr, spacings_f, Dparallel_f)#, maxfev=thismaxfev)
sigma_rm_f           = popt[0]
# Ds
Ds_rm_rw   = catex_resin_membr(spacings_rw,sigma_rm_rw)
Ds_rm_stat = catex_resin_membr(spacings_stat,sigma_rm_stat)
Ds_rm_dyn  = catex_resin_membr(spacings_dyn,sigma_rm_dyn)
Ds_rm_f    = catex_resin_membr(spacings_f,sigma_rm_f)

## Plotting:# Need fancy plotting. Find code and redo fancy plotting.
plt.figure(figsize=(17,5))
ax = plt.subplot(111)
ax.errorbar(spacings_dyn, Dparallel_dyn, yerr=Dparallel_stdv_dyn, color='g', capsize=2, label=r'$D_\parallel$, dyn.')
#ax.errorbar(spacings_stat, Dparallel_stat, yerr=Dparallel_stdv_stat, color='limegreen', capsize=2, label=r'$D_\parallel$, stat.')
#ax.errorbar(spacings_f, Dparallel_f, yerr=Dparallel_stdv_f, color='springgreen', capsize=2, label=r'$D_\parallel$, straight')
#ax.plot(dg,DRW, color='b', label='Random walk, fixed geom., no refl.')
ax.plot(spacings_dyn,Ds_op_dyn, '--', label=r'Ordered packings, $\sigma=$%.2f' % sigma_op_dyn) # I have most values for the dynamic brush
ax.plot(spacings_dyn,Ds_hr_dyn, '--', label=r'Hyperbola of revolution, $\sigma=$%.2f' % sigma_hr_dyn) # 2
#ax.plot(spacings_dyn,Ds_nm_dyn, '--', label=r'Not monosized spheres, $\sigma=$%.2f' % sigma_nm_dyn) 
#ax.plot(spacings_dyn,Ds_ps_dyn, '--', label=r'Monodisperse sphere packing, $\sigma=$%.2f' % sigma_ps_dyn)
#ax.plot(spacings_dyn,Ds_os_dyn, '--', label=r'Overlapping spheres, $\sigma=$%.2f' % sigma_os_dyn) # 5
#ax.plot(spacings_dyn,Ds_oc_dyn, '--', label=r'Overlapping cylinders, $\sigma=$%.2f' % sigma_oc_dyn) 
ax.plot(spacings_dyn,Ds_hc_dyn, '--', label=r'Heterogeneous catalyst, $\sigma=$%.2f' % sigma_hc_dyn) 
ax.plot(spacings_dyn,Ds_rm_dyn, '--', label=r'Cation-exchange resin memrane, $\sigma=$%.2f' % sigma_rm_dyn) 
# Add more plotting! # Need to beware of negative porosity for some!
plt.xlabel(r'Spacing $d$')
plt.ylabel(r'Normalized diffusion coefficient $D/D_{bulk}$')
plt.title(r'$D/D_{bulk}$ vs $d$')
#plt.legend(loc='lower right', prop={'size': 12})
#plt.tight_layout()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_d_dyn)


## Plotting:# Need fancy plotting. Find code and redo fancy plotting.
plt.figure(figsize=(17,5))
ax = plt.subplot(111)
#ax.errorbar(spacings_dyn, Dparallel_dyn, yerr=Dparallel_stdv_dyn, color='g', capsize=2, label=r'$D_\parallel$, dyn.')
ax.errorbar(spacings_stat, Dparallel_stat, yerr=Dparallel_stdv_stat, color='limegreen', capsize=2, label=r'$D_\parallel$, stat.')
#ax.errorbar(spacings_f, Dparallel_f, yerr=Dparallel_stdv_f, color='springgreen', capsize=2, label=r'$D_\parallel$, straight')
#ax.plot(dg,DRW, color='b', label='Random walk, fixed geom., no refl.')
ax.plot(spacings_stat,Ds_op_stat, '--', label=r'Ordered packings, $\sigma=$%.2f' % sigma_op_stat) # I have most values for the dynamic brush
ax.plot(spacings_stat,Ds_hr_stat, '--', label=r'Hyperbola of revolution, $\sigma=$%.2f' % sigma_hr_stat) # 2
#ax.plot(spacings_stat,Ds_nm_stat, '--', label=r'Not monosized spheres, $\sigma=$%.2f' % sigma_nm_stat) 
#ax.plot(spacings_stat,Ds_ps_stat, '--', label=r'Monodisperse sphere packing, $\sigma=$%.2f' % sigma_ps_stat)
#ax.plot(spacings_stat,Ds_os_stat, '--', label=r'Overlapping spheres, $\sigma=$%.2f' % sigma_os_stat) # 5
#ax.plot(spacings_stat,Ds_oc_stat, '--', label=r'Overlapping cylinders, $\sigma=$%.2f' % sigma_oc_stat) 
ax.plot(spacings_stat,Ds_hc_stat, '--', label=r'Heterogeneous catalyst, $\sigma=$%.2f' % sigma_hc_stat) 
ax.plot(spacings_stat,Ds_rm_stat, '--', label=r'Cation-exchange resin memrane, $\sigma=$%.2f' % sigma_rm_stat) 
# Add more plotting! # Need to beware of negative porosity for some!
plt.xlabel(r'Spacing $d$')
plt.ylabel(r'Normalized diffusion coefficient $D/D_{bulk}$')
plt.title(r'$D/D_{bulk}$ vs $d$')
#plt.legend(loc='lower right', prop={'size': 12})
#plt.tight_layout()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_d_stat)




## Plotting:# Need fancy plotting. Find code and redo fancy plotting.
plt.figure(figsize=(17,5))
ax = plt.subplot(111)
#ax.errorbar(spacings_dyn, Dparallel_dyn, yerr=Dparallel_stdv_dyn, color='g', capsize=2, label=r'$D_\parallel$, dyn.')
#ax.errorbar(spacings_stat, Dparallel_stat, yerr=Dparallel_stdv_stat, color='limegreen', capsize=2, label=r'$D_\parallel$, stat.')
ax.errorbar(spacings_f, Dparallel_f, yerr=Dparallel_stdv_f, color='springgreen', capsize=2, label=r'$D_\parallel$, straight')
#ax.plot(dg,DRW, color='b', label='Random walk, fixed geom., no refl.')
ax.plot(spacings_f,Ds_op_f, '--', label=r'Ordered packings, $\sigma=$%.2f' % sigma_op_f) # I have most values for the dynamic brush
ax.plot(spacings_f,Ds_hr_f, '--', label=r'Hyperbola of revolution, $\sigma=$%.2f' % sigma_hr_f) # 2
#ax.plot(spacings_f,Ds_nm_f, '--', label=r'Not monosized spheres, $\sigma=$%.2f' % sigma_nm_f) 
#ax.plot(spacings_f,Ds_ps_f, '--', label=r'Monodisperse sphere packing, $\sigma=$%.2f' % sigma_ps_f)
#ax.plot(spacings_f,Ds_os_f, '--', label=r'Overlapping spheres, $\sigma=$%.2f' % sigma_os_f) # 5
#ax.plot(spacings_f,Ds_oc_f, '--', label=r'Overlapping cylinders, $\sigma=$%.2f' % sigma_oc_f) 
ax.plot(spacings_f,Ds_hc_f, '--', label=r'Heterogeneous catalyst, $\sigma=$%.2f' % sigma_hc_f) 
ax.plot(spacings_f,Ds_rm_f, '--', label=r'Cation-exchange resin memrane, $\sigma=$%.2f' % sigma_rm_f) 
# Add more plotting! # Need to beware of negative porosity for some!
plt.xlabel(r'Spacing $d$')
plt.ylabel(r'Normalized diffusion coefficient $D/D_{bulk}$')
plt.title(r'$D/D_{bulk}$ vs $d$')
#plt.legend(loc='lower right', prop={'size': 12})
#plt.tight_layout()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_d_f)

## RMSDs:
# Ordered packing: 1
D_rmsd_op_rw   = rmsd(DRW,Ds_op_rw)
D_rmsd_op_stat = rmsd(Dparallel_stat,Ds_op_stat)
D_rmsd_op_dyn  = rmsd(Dparallel_dyn,Ds_op_dyn)
D_rmsd_op_f    = rmsd(Dparallel_f,Ds_op_f)
# Hyperbola of revolution 2
D_rmsd_hr_rw   = rmsd(DRW,Ds_hr_rw)
D_rmsd_hr_stat = rmsd(Dparallel_stat,Ds_hr_stat)
D_rmsd_hr_dyn  = rmsd(Dparallel_dyn,Ds_hr_dyn)
D_rmsd_hr_f    = rmsd(Dparallel_f,Ds_hr_f)
# Not monosized spheres 3
D_rmsd_nm_rw   = rmsd(DRW,Ds_nm_rw)      # Used to need _D_rw_positive
D_rmsd_nm_stat = rmsd(Dparallel_stat,Ds_nm_stat)
D_rmsd_nm_dyn  = rmsd(Dparallel_dyn,Ds_nm_dyn)
D_rmsd_nm_f    = rmsd(Dparallel_f,Ds_nm_f)
# pshimd_spherepacking 4
D_rmsd_ps_rw   = rmsd(DRW,Ds_ps_rw)      # Used to need _positive
D_rmsd_ps_stat = rmsd(Dparallel_stat,Ds_ps_stat)
D_rmsd_ps_dyn  = rmsd(Dparallel_dyn,Ds_ps_dyn)
D_rmsd_ps_f    = rmsd(Dparallel_f,Ds_ps_f)
# Overlapping spheres 5
D_rmsd_os_rw   = rmsd(DRW,Ds_os_rw)      # Used to need _D_rw_positive
D_rmsd_os_stat = rmsd(Dparallel_stat,Ds_os_stat)
D_rmsd_os_dyn  = rmsd(Dparallel_dyn,Ds_os_dyn)
D_rmsd_os_f    = rmsd(Dparallel_f,Ds_os_f)
# Overlapping cylinders 6
D_rmsd_oc_rw   = rmsd(DRW,Ds_oc_rw)      # Used to need _D_rw_positive
D_rmsd_oc_stat = rmsd(Dparallel_stat,Ds_oc_stat)
D_rmsd_oc_dyn  = rmsd(Dparallel_dyn,Ds_oc_dyn)
D_rmsd_oc_f    = rmsd(Dparallel_f,Ds_oc_f)
# Heterogeneous catalyst 7
D_rmsd_hc_rw   = rmsd(DRW,Ds_hc_rw)
D_rmsd_hc_stat = rmsd(Dparallel_stat,Ds_hc_stat)
D_rmsd_hc_dyn  = rmsd(Dparallel_dyn,Ds_hc_dyn)
D_rmsd_hc_f    = rmsd(Dparallel_f,Ds_hc_f)
# Cation exchange resin memrane 8
D_rmsd_rm_rw   = rmsd(DRW,Ds_rm_rw)
D_rmsd_rm_stat = rmsd(Dparallel_stat,Ds_rm_stat)
D_rmsd_rm_dyn  = rmsd(Dparallel_dyn,Ds_rm_dyn)
D_rmsd_rm_f    = rmsd(Dparallel_f,Ds_rm_f)


# Write rmsd to file
'''
rmsdfilename_RW
rmsdfilename_forest
rmsdfilename_dyn
rmsdfilename_stat
'''

rmsdfile_rw   = open(rmsdfilename_RW,'w')
rmsdfile_f    = open(rmsdfilename_forest,'w')
rmsdfile_dyn  = open(rmsdfilename_dyn,'w')
rmsdfile_stat = open(rmsdfilename_stat,'w')
# Write a header!!! # Or have rows: name value? # One line is easier to read in.
#                                                                        1           2            3              4          5             6           7             8
rmsdfile_rw.write('%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n' % (D_rmsd_op_rw,D_rmsd_hr_rw,D_rmsd_nm_rw,D_rmsd_ps_rw,D_rmsd_os_rw,D_rmsd_oc_rw,D_rmsd_hc_rw,D_rmsd_rm_rw))
rmsdfile_rw.write('sigma: %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n' % (sigma_op_rw,sigma_hr_rw,sigma_nm_rw,sigma_ps_rw,sigma_os_rw,sigma_oc_rw,sigma_hc_rw,sigma_rm_rw))
rmsdfile_rw.close()

rmsdfile_f = open(rmsdfilename_forest,'w')
rmsdfile_f.write('%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n' % (D_rmsd_op_f,D_rmsd_hr_f,D_rmsd_nm_f,D_rmsd_ps_f,D_rmsd_os_f,D_rmsd_oc_f,D_rmsd_hc_f,D_rmsd_rm_f))
rmsdfile_f.write('sigma: %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n' % (sigma_op_f,sigma_hr_f,sigma_nm_f,sigma_ps_f,sigma_os_f,sigma_oc_f,sigma_hc_f,sigma_rm_f))
rmsdfile_f.close()

rmsdfile_dyn = open(rmsdfilename_dyn,'w')
rmsdfile_dyn.write('%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n' % (D_rmsd_op_dyn,D_rmsd_hr_dyn,D_rmsd_nm_dyn,D_rmsd_ps_dyn,D_rmsd_os_dyn,D_rmsd_oc_dyn,D_rmsd_hc_dyn,D_rmsd_rm_dyn))
rmsdfile_dyn.write('sigma: %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n' % (sigma_op_dyn,sigma_hr_dyn,sigma_nm_dyn,sigma_ps_dyn,sigma_os_dyn,sigma_oc_dyn,sigma_hc_dyn,sigma_rm_dyn))
rmsdfile_dyn.close()

rmsdfile_stat = open(rmsdfilename_stat,'w')
rmsdfile_stat.write('%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n' % (D_rmsd_op_stat,D_rmsd_hr_stat,D_rmsd_nm_stat,D_rmsd_ps_stat,D_rmsd_os_stat,D_rmsd_oc_stat,D_rmsd_hc_stat,D_rmsd_rm_stat))
rmsdfile_stat.write('sigma: %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n' % (sigma_op_stat,sigma_hr_stat,sigma_nm_stat,sigma_ps_stat,sigma_os_stat,sigma_oc_stat,sigma_hc_stat,sigma_rm_stat))
rmsdfile_stat.close()

if showplots==True:
    plt.show()
