import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob

thismaxfev = 10000

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

def mymodel1(d,k,sigma):
    phi = 1 - np.pi*(sigma/d)
    teller = 1-phi*k
    nevner = phi*(1-k)
    tau = (teller/nevner)**2
    return 1./tau

def mymodel2(d,k,sigma):
    phi = 1 - np.pi*(sigma/d)
    teller = 1+(1-phi)*(1-k)
    nevner = 1-(1-phi)*k
    tau = (teller/nevner)**2
    return 1./tau

def mymodel3(d,k,sigma):
    phi = 1 - np.pi*(sigma/d)
    teller = 1-phi
    nevner = 1-k
    tau = (1+teller/nevner)**2
    return 1./tau

def mymodel4(d,k,f,sigma):
    # Assumes probability of neighbor being occupied when you're occupied: kvp**f
    phi = 1 - np.pi*(sigma/d)
    vp = 1-phi
    nevner = 1-k*vp**f
    tau = (1+vp/nevner)**2
    return 1./tau

def mymodel5(d,k,f,sigma):
    # Assumes probability of neighbor being occupied when you're occupied: kvp**f
    phi = 1 - np.pi*(sigma/d)
    vp = 1-phi
    tau = (1+vp/(1-(k+((k-1)*vp**f))))**2
    return 1./tau

def mymodel6(d,k,sigma):
    # Assumes probability of neighbor being occupied when you're occupied: kvp**f
    phi = 1 - np.pi*(sigma/d) 
    vp = 1-phi
    nevner = 1-k*vp**d
    tau = (1+vp/nevner)**2
    return 1./tau

def mymodel7(d,k,sigma):
    # Assumes probability of neighbor being occupied when you're occupied: kvp**f
    phi = 1 - np.pi*(sigma/d) 
    vp = 1-phi
    nevner = 1-k*vp**(-1./d)
    tau = (1+vp/nevner)**2
    return 1./tau

def powerlaw1(d,k,l,n):
    return 1 - (d/l)**n*k

# Misc.
def giveporosity(d,sigma):
    phi = 1 - np.pi*(sigma/d)
    #print('phi:', phi, ' when sigma=',sigma)
    return phi

def findfirstpositive(x):
    haspassed = False
    N = len(x)
    for i in range(N):
        if haspassed==False and x[i]>0:
            haspassed=True
            index=i
    if haspassed==False:
        i=N
    return i

def find_d_for_first_positive(x,d):
    haspassed = False
    N = len(x)
    for i in range(N):
        if haspassed==False and x[i]>0:
            haspassed=True
            index=i
            passd = d[i]
    if haspassed==False:
        passd='--'    
    return passd

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
fitds         = np.linspace(2,25,100)

endlocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/Nocut/'
if moresigmas==True:
    endlocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Moresigmas/Nocut/'

basepath_base      = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/'
endlocation_static = basepath_base+'D_vs_d/Nocut/'

bulklocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/'

bulkfilename  = bulklocation + 'diffusion_bulk'+filestext
if bulk_cut==True:
    bulkfilename = bulkfilename +'_cut.txt'
else:
    bulkfilename = bulkfilename +'_uncut.txt'

## Files to read
brushfilename_dyn  = endlocation + 'D_vs_d_better_rms_Nestimates10.txt'
brushfilename_stat = endlocation_static + 'D_vs_d_static_better_rms_Nestimates10.txt'

## Files to write to
rmsdfilename_RW = endlocation_static+'D_rmsd_close_packing_vs_RWgeom_norefl_rsphere'+str(rsphere)+'_maxlen1.0_modelr_findsigma_fittomid.txt'
rmsdfilename_forest = endlocation_static+'D_rmsd_close_packing_vs_forest_modelr_findsigma_fittomid.txt'
rmsdfilename_dyn = endlocation_static+'D_rmsd_close_packing_vs_dyn_modelr_findsigma_fittomid.txt'
rmsdfilename_stat = endlocation_static+'D_rmsd_close_packing_vs_stat_modelr_findsigma_fittomid.txt'
###
if big==False:
    ### d
    plotname_d_dyn     = endlocation_static+'D_vs_d_dyn_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_fittomid_sigmafit_upperbounds1inf.png'
    plotname_d_stat    = endlocation_static+'D_vs_d_stat_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_fittomid_sigmafit_upperbounds1inf.png'
    plotname_d_f       = endlocation_static+'D_vs_d_f_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_fittomid_sigmafit_upperbounds1inf.png'
    # perp:
    plotname_d_dyn_perp     = endlocation_static+'D_vs_d_dyn_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_perp_fittomid_sigmafit_upperbounds1inf.png'
    plotname_d_stat_perp    = endlocation_static+'D_vs_d_stat_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_perp_fittomid_sigmafit_upperbounds1inf.png'
    plotname_d_f_perp       = endlocation_static+'D_vs_d_f_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_perp_fittomid_sigmafit_upperbounds1inf.png'
    plotname_cut_d = endlocation_static+'D_vs_d_dyn_vs_stat_vs_RWgeom_cut_norefl_rsphere'+str(rsphere)+'_modelr_findsigma_fittomid_sigmafit_upperbounds1inf.png' ## Use this?
else:
    ### d 
    plotname_d_dyn     = endlocation_static+'D_vs_d_dyn_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_fittomid_sigmafit_upperbounds1inf.png'
    plotname_d_stat    = endlocation_static+'D_vs_d_stat_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_fittomid_sigmafit_upperbounds1inf.png'
    plotname_d_f       = endlocation_static+'D_vs_d_f_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_fittomid_sigmafit_upperbounds1inf.png'
    # perp:
    plotname_d_dyn_perp     = endlocation_static+'D_vs_d_dyn_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_perp_fittomid_sigmafit_upperbounds1inf.png'
    plotname_d_stat_perp    = endlocation_static+'D_vs_d_stat_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_perp_fittomid_sigmafit_upperbounds1inf.png'
    plotname_d_f_perp       = endlocation_static+'D_vs_d_f_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_perp_fittomid_sigmafit_upperbounds1inf.png'
    plotname_cut_d = endlocation_static+'D_vs_d_dyn_vs_stat_cut_vs_RWgeom_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_fittomid_sigmafit_upperbounds1inf.png' ## Use this?

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
endlocation_forest = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_forest/D_vs_d/Nocut/'
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
folder   = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Randomwalk_fixedgeometry/'

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


indices_dyn = []
for i in range(N_dyn):
    if spacings_dyn[i]>=2 and spacings_dyn[i]<30:
        indices_dyn.append(i)
startindex_dyn = indices_dyn[0]
endindex_dyn   = indices_dyn[-1]+1
spacings_dyn   = spacings_dyn[startindex_dyn:endindex_dyn]
Dparallel_dyn  = Dparallel_dyn[startindex_dyn:endindex_dyn]
Dparallel_stdv_dyn  = Dparallel_stdv_dyn[startindex_dyn:endindex_dyn]
Dzs_dyn        = Dzs_dyn[startindex_dyn:endindex_dyn]
Dzs_stdv_dyn   = Dzs_stdv_dyn[startindex_dyn:endindex_dyn]

indices_stat = []
for i in range(N_stat):
    if spacings_stat[i]>=2 and spacings_stat[i]<30:
        indices_stat.append(i)
startindex_stat = indices_stat[0]
endindex_stat   = indices_stat[-1]+1
spacings_stat   = spacings_stat[startindex_stat:endindex_stat]
Dparallel_stat  = Dparallel_stat[startindex_stat:endindex_stat]
Dparallel_stdv_stat  = Dparallel_stdv_stat[startindex_stat:endindex_stat]
Dzs_stat        = Dzs_stat[startindex_stat:endindex_stat]
Dzs_stdv_stat   = Dzs_stdv_stat[startindex_stat:endindex_stat]

indices_f = []
for i in range(N_f):
    if spacings_f[i]>=2 and spacings_f[i]<30:
        indices_f.append(i)
startindex_f = indices_f[0]
endindex_f   = indices_f[-1]+1
spacings_f   = spacings_f[startindex_f:endindex_f]
Dparallel_f  = Dparallel_f[startindex_f:endindex_f]
Dparallel_stdv_f  = Dparallel_stdv_f[startindex_f:endindex_f]
Dzs_f  = Dzs_f[startindex_f:endindex_f]
Dzs_stdv_f  = Dzs_stdv_f[startindex_f:endindex_f]

# Porosity
dg = np.array(dg)
spacings_rw   = dg
spacings_stat = np.array(spacings_stat)
spacings_dyn  = np.array(spacings_dyn)
spacings_f    = np.array(spacings_f)

# Ordered packing
popt, pcov           = curve_fit(ordered_packings, spacings_rw, DRW, maxfev=thismaxfev)
sigma_op_rw          = popt[0]
popt, pcov           = curve_fit(ordered_packings, spacings_stat, Dparallel_stat, maxfev=thismaxfev)
sigma_op_stat        = popt[0]
sigma_op_stat_rms    = np.sqrt(np.diag(pcov))[0]
popt, pcov           = curve_fit(ordered_packings, spacings_dyn, Dparallel_dyn, p0=7, maxfev=thismaxfev)
sigma_op_dyn         = popt[0]
sigma_op_dyn_rms     = np.sqrt(np.diag(pcov))[0]
popt, pcov           = curve_fit(ordered_packings, spacings_f, Dparallel_f, maxfev=thismaxfev)
sigma_op_f           = popt[0]
# Perp
popt, pcov             = curve_fit(ordered_packings, spacings_stat, Dzs_stat, maxfev=thismaxfev)
sigma_op_stat_perp     = popt[0]
sigma_op_stat_perp_rms = np.sqrt(np.diag(pcov))[0]
popt, pcov             = curve_fit(ordered_packings, spacings_dyn, Dzs_dyn, p0=7, maxfev=thismaxfev)
sigma_op_dyn_perp      = popt[0]
sigma_op_dyn_perp_rms  = np.sqrt(np.diag(pcov))[0]
popt, pcov             = curve_fit(ordered_packings, spacings_f, Dzs_f, maxfev=thismaxfev)
sigma_op_f_perp        = popt[0]

#sigma_op_dyn = 7
# Ds
Ds_op_rw   = ordered_packings(fitds,sigma_op_rw)
Ds_op_stat = ordered_packings(fitds,sigma_op_stat)
Ds_op_dyn  = ordered_packings(fitds,sigma_op_dyn)
Ds_op_f    = ordered_packings(fitds,sigma_op_f)
Ds_op_stat_perp = ordered_packings(fitds,sigma_op_stat_perp)
Ds_op_dyn_perp  = ordered_packings(fitds,sigma_op_dyn_perp)
Ds_op_f_perp    = ordered_packings(fitds,sigma_op_f_perp)

# Hyperbola of revolution
popt, pcov           = curve_fit(hyp_rev, spacings_rw, DRW, maxfev=thismaxfev)
sigma_hr_rw          = popt[0]
popt, pcov           = curve_fit(hyp_rev, spacings_stat, Dparallel_stat, maxfev=thismaxfev)
sigma_hr_stat        = popt[0]
sigma_hr_stat_rms    = np.sqrt(np.diag(pcov))[0]
popt, pcov           = curve_fit(hyp_rev, spacings_dyn, Dparallel_dyn, p0=5, maxfev=thismaxfev)
sigma_hr_dyn         = popt[0]
sigma_hr_dyn_rms     = np.sqrt(np.diag(pcov))[0]
popt, pcov           = curve_fit(hyp_rev, spacings_f, Dparallel_f, maxfev=thismaxfev)
sigma_hr_f           = popt[0]
# Perp
popt, pcov             = curve_fit(hyp_rev, spacings_stat, Dzs_stat, maxfev=thismaxfev)
sigma_hr_stat_perp     = popt[0]
sigma_hr_stat_perp_rms = np.sqrt(np.diag(pcov))[0]
popt, pcov             = curve_fit(hyp_rev, spacings_dyn, Dzs_dyn, p0=5, maxfev=thismaxfev)
sigma_hr_dyn_perp      = popt[0]
sigma_hr_dyn_perp_rms  = np.sqrt(np.diag(pcov))[0]
popt, pcov             = curve_fit(hyp_rev, spacings_f, Dzs_f, maxfev=thismaxfev)
sigma_hr_f_perp        = popt[0]

#sigma_hr_dyn = 5
# Ds
Ds_hr_rw   = hyp_rev(fitds,sigma_hr_rw)
Ds_hr_stat = hyp_rev(fitds,sigma_hr_stat)
Ds_hr_dyn  = hyp_rev(fitds,sigma_hr_dyn)
Ds_hr_f    = hyp_rev(fitds,sigma_hr_f)
Ds_hr_stat_perp = hyp_rev(fitds,sigma_hr_stat_perp)
Ds_hr_dyn_perp  = hyp_rev(fitds,sigma_hr_dyn_perp)
Ds_hr_f_perp    = hyp_rev(fitds,sigma_hr_f_perp)

# Not monosized spheres
'''
popt, pcov           = curve_fit(notmonosph, spacings_rw, DRW, maxfev=thismaxfev)
sigma_nm_rw          = popt[0]
popt, pcov           = curve_fit(notmonosph, spacings_stat, Dparallel_stat, p0=0.5, maxfev=thismaxfev)
sigma_nm_stat        = popt[0]
popt, pcov           = curve_fit(notmonosph, spacings_dyn, Dparallel_dyn, p0=0.5, maxfev=thismaxfev)
sigma_nm_dyn         = popt[0]
popt, pcov           = curve_fit(notmonosph, spacings_f, Dparallel_f, p0=0.5, maxfev=thismaxfev)
sigma_nm_f           = popt[0]
# Ds
Ds_nm_rw   = notmonosph(fitds,sigma_nm_rw)
Ds_nm_stat = notmonosph(fitds,sigma_nm_stat)
Ds_nm_dyn  = notmonosph(fitds,sigma_nm_dyn)
Ds_nm_f    = notmonosph(fitds,sigma_nm_f)
# pshimd_spherepacking
popt, pcov           = curve_fit(pshimd_spherepacking, spacings_rw, DRW, maxfev=thismaxfev)
sigma_ps_rw          = popt[0]
popt, pcov           = curve_fit(pshimd_spherepacking, spacings_stat, Dparallel_stat, p0=0.5, maxfev=thismaxfev)
sigma_ps_stat        = popt[0]
popt, pcov           = curve_fit(pshimd_spherepacking, spacings_dyn, Dparallel_dyn, p0=0.5, maxfev=thismaxfev)
sigma_ps_dyn         = popt[0]
popt, pcov           = curve_fit(pshimd_spherepacking, spacings_f, Dparallel_f, p0=0.5, maxfev=thismaxfev)
sigma_ps_f           = popt[0]

# Ds
Ds_ps_rw   = pshimd_spherepacking(fitds,sigma_ps_rw)
Ds_ps_stat = pshimd_spherepacking(fitds,sigma_ps_stat)
Ds_ps_dyn  = pshimd_spherepacking(fitds,sigma_ps_dyn)
Ds_ps_f    = pshimd_spherepacking(fitds,sigma_ps_f)
# Overlapping spheres
popt, pcov           = curve_fit(overlapping_spheres, spacings_rw, DRW, maxfev=thismaxfev)
sigma_os_rw          = popt[0]
popt, pcov           = curve_fit(overlapping_spheres, spacings_stat, Dparallel_stat, p0=0.5, maxfev=thismaxfev)
sigma_os_stat        = popt[0]
popt, pcov           = curve_fit(overlapping_spheres, spacings_dyn, Dparallel_dyn,p0=0.5, maxfev=thismaxfev)
sigma_os_dyn         = popt[0]
popt, pcov           = curve_fit(overlapping_spheres, spacings_f, Dparallel_f, p0=0.5, maxfev=thismaxfev)
sigma_os_f           = popt[0]
# Ds
Ds_os_rw   = overlapping_spheres(fitds,sigma_os_rw)
Ds_os_stat = overlapping_spheres(fitds,sigma_os_stat)
Ds_os_dyn  = overlapping_spheres(fitds,sigma_os_dyn)
Ds_os_f    = overlapping_spheres(fitds,sigma_os_f)
# Overlapping cylinders
popt, pcov           = curve_fit(overlapping_cylinders, spacings_rw, DRW, maxfev=thismaxfev)
sigma_oc_rw          = popt[0]
popt, pcov           = curve_fit(overlapping_cylinders, spacings_stat, Dparallel_stat, sigma=sigmas_stat, p0=0.5, maxfev=thismaxfev)
sigma_oc_stat        = popt[0]
popt, pcov           = curve_fit(overlapping_cylinders, spacings_dyn, Dparallel_dyn, p0=0.5, maxfev=thismaxfev)
sigma_oc_dyn         = popt[0]
popt, pcov           = curve_fit(overlapping_cylinders, spacings_f, Dparallel_f, p0=0.5, maxfev=thismaxfev)
sigma_oc_f           = popt[0]
# Ds
Ds_oc_rw   = overlapping_cylinders(fitds,sigma_oc_rw)
Ds_oc_stat = overlapping_cylinders(fitds,sigma_oc_stat)
Ds_oc_dyn  = overlapping_cylinders(fitds,sigma_oc_dyn)
Ds_oc_f    = overlapping_cylinders(fitds,sigma_oc_f)
'''

# Heterogeneous catalyst
popt, pcov           = curve_fit(het_cat, spacings_rw, DRW, maxfev=thismaxfev)
sigma_hc_rw          = popt[0]
popt, pcov           = curve_fit(het_cat, spacings_stat, Dparallel_stat,maxfev=thismaxfev)
sigma_hc_stat        = popt[0]
sigma_hc_stat_rms    = np.sqrt(np.diag(pcov))[0]
popt, pcov           = curve_fit(het_cat, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev)
sigma_hc_dyn         = popt[0]
sigma_hc_dyn_rms     = np.sqrt(np.diag(pcov))[0]
popt, pcov           = curve_fit(het_cat, spacings_f, Dparallel_f, maxfev=thismaxfev)
sigma_hc_f           = popt[0]
# Perp
popt, pcov             = curve_fit(het_cat, spacings_stat, Dzs_stat,maxfev=thismaxfev)
sigma_hc_stat_perp     = popt[0]
sigma_hc_stat_perp_rms = np.sqrt(np.diag(pcov))[0]
popt, pcov             = curve_fit(het_cat, spacings_dyn, Dzs_dyn, maxfev=thismaxfev)
sigma_hc_dyn_perp      = popt[0]
sigma_hc_dyn_perp_rms  = np.sqrt(np.diag(pcov))[0]
popt, pcov             = curve_fit(het_cat, spacings_f, Dzs_f, maxfev=thismaxfev)
sigma_hc_f_perp        = popt[0]

# Ds
Ds_hc_rw   = het_cat(fitds,sigma_hc_rw)
Ds_hc_stat = het_cat(fitds,sigma_hc_stat)
Ds_hc_dyn  = het_cat(fitds,sigma_hc_dyn)
Ds_hc_f    = het_cat(fitds,sigma_hc_f)
Ds_hc_stat_perp = het_cat(fitds,sigma_hc_stat_perp)
Ds_hc_dyn_perp  = het_cat(fitds,sigma_hc_dyn_perp)
Ds_hc_f_perp    = het_cat(fitds,sigma_hc_f_perp)

# Cation exchange resin membrane
popt, pcov           = curve_fit(catex_resin_membr, spacings_rw, DRW, maxfev=thismaxfev)
sigma_rm_rw          = popt[0]
popt, pcov           = curve_fit(catex_resin_membr, spacings_stat, Dparallel_stat, maxfev=thismaxfev)
sigma_rm_stat        = popt[0]
sigma_rm_stat_rms    = np.sqrt(np.diag(pcov))[0]
popt, pcov           = curve_fit(catex_resin_membr, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev)
sigma_rm_dyn         = popt[0]
sigma_rm_dyn_rms     = np.sqrt(np.diag(pcov))[0]
popt, pcov           = curve_fit(catex_resin_membr, spacings_f, Dparallel_f, maxfev=thismaxfev)
sigma_rm_f           = popt[0]
#Perp
popt, pcov             = curve_fit(catex_resin_membr, spacings_stat, Dzs_stat, maxfev=thismaxfev)
sigma_rm_stat_perp     = popt[0]
sigma_rm_stat_perp_rms = np.sqrt(np.diag(pcov))[0]
popt, pcov             = curve_fit(catex_resin_membr, spacings_dyn, Dzs_dyn, maxfev=thismaxfev)
sigma_rm_dyn_perp      = popt[0]
sigma_rm_dyn_perp_rms  = np.sqrt(np.diag(pcov))[0]
popt, pcov             = curve_fit(catex_resin_membr, spacings_f, Dzs_f, maxfev=thismaxfev)
sigma_rm_f_perp        = popt[0]

# Ds
Ds_rm_rw   = catex_resin_membr(fitds,sigma_rm_rw)
Ds_rm_stat = catex_resin_membr(fitds,sigma_rm_stat)
Ds_rm_dyn  = catex_resin_membr(fitds,sigma_rm_dyn)
Ds_rm_f    = catex_resin_membr(fitds,sigma_rm_f)
Ds_rm_stat_perp = catex_resin_membr(fitds,sigma_rm_stat_perp)
Ds_rm_dyn_perp  = catex_resin_membr(fitds,sigma_rm_dyn_perp)
Ds_rm_f_perp    = catex_resin_membr(fitds,sigma_rm_f_perp)

# My models:
### mymodel1(d,k,sigma)
popt, pcov         = curve_fit(mymodel1, spacings_rw, DRW, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
k_m1_rw            = popt[0]
sigma_m1_rw        = popt[1]
popt, pcov         = curve_fit(mymodel1, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m1_stat          = popt[0]
sigma_m1_stat      = popt[1]
k_m1_stat_stdv     = rms_params[0]
sigma_m1_stat_stdv = rms_params[1]
popt, pcov         = curve_fit(mymodel1, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m1_dyn           = popt[0]
sigma_m1_dyn       = popt[1]
k_m1_dyn_stdv      = rms_params[0]
sigma_m1_dyn_stdv  = rms_params[1]
popt, pcov         = curve_fit(mymodel1, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m1_f             = popt[0]
sigma_m1_f         = popt[1]
k_m1_dyn_stdv      = rms_params[0]
sigma_m1_dyn_stdv  = rms_params[1]
#Perp
popt, pcov              = curve_fit(mymodel1, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m1_stat_perp          = popt[0]
k_m1_stat_perp_stdv     = rms_params[0]
sigma_m1_stat_perp      = popt[1]
sigma_m1_stat_perp_stdv = rms_params[1]
popt, pcov              = curve_fit(mymodel1, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m1_dyn_perp           = popt[0]
k_m1_dyn_perp_stdv      = rms_params[0]
sigma_m1_dyn_perp       = popt[1]
sigma_m1_dyn_perp_stdv  = rms_params[1]
popt, pcov              = curve_fit(mymodel1, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m1_f_perp             = popt[0]
sigma_m1_f_perp         = popt[1]

# Ds
Ds_m1_rw   = mymodel1(fitds,k_m1_rw,sigma_m1_rw)
Ds_m1_stat = mymodel1(fitds,k_m1_stat,sigma_m1_stat)
Ds_m1_dyn  = mymodel1(fitds,k_m1_dyn,sigma_m1_dyn)
Ds_m1_f    = mymodel1(fitds,k_m1_f,sigma_m1_f)
Ds_m1_stat_perp = mymodel1(fitds,k_m1_stat_perp,sigma_m1_stat_perp)
Ds_m1_dyn_perp  = mymodel1(fitds,k_m1_dyn_perp,sigma_m1_dyn_perp)
Ds_m1_f_perp    = mymodel1(fitds,k_m1_f_perp,sigma_m1_f_perp)

### mymodel2(d,k)
popt, pcov         = curve_fit(mymodel2, spacings_rw, DRW, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
k_m2_rw            = popt[0]
sigma_m2_rw        = popt[1]
popt, pcov = curve_fit(mymodel2, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m2_stat          = popt[0]
k_m2_stat_stdv     = rms_params[0]
sigma_m2_stat      = popt[1]
sigma_m2_stat_stdv = rms_params[1]
popt, pcov         = curve_fit(mymodel2, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m2_dyn           = popt[0]
k_m2_dyn_stdv      = rms_params[0]
sigma_m2_dyn       = popt[1]
sigma_m2_dyn_stdv  = rms_params[1]
popt, pcov         = curve_fit(mymodel2, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
k_m2_f             = popt[0]
sigma_m2_f         = popt[1]
#Perp
popt, pcov   = curve_fit(mymodel2, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m2_stat_perp          = popt[0]
k_m2_stat_perp_stdv     = rms_params[0]
sigma_m2_stat_perp      = popt[1]
sigma_m2_stat_perp_stdv = rms_params[1]
popt, pcov              = curve_fit(mymodel2, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m2_dyn_perp           = popt[0]
k_m2_dyn_perp_stdv      = rms_params[0]
sigma_m2_dyn_perp       = popt[1]
sigma_m2_dyn_perp_stdv  = rms_params[1]
popt, pcov              = curve_fit(mymodel2, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
k_m2_f_perp             = popt[0]
sigma_m2_f_perp         = popt[1]

# Ds
Ds_m2_rw   = mymodel2(fitds,k_m2_rw,sigma_m2_rw)
Ds_m2_stat = mymodel2(fitds,k_m2_stat,sigma_m2_stat)
Ds_m2_dyn  = mymodel2(fitds,k_m2_dyn,sigma_m2_dyn)
Ds_m2_f    = mymodel2(fitds,k_m2_f,sigma_m2_f)
Ds_m2_stat_perp = mymodel2(fitds,k_m2_stat_perp,sigma_m2_stat_perp)
Ds_m2_dyn_perp  = mymodel2(fitds,k_m2_dyn_perp,sigma_m2_dyn_perp)
Ds_m2_f_perp    = mymodel2(fitds,k_m2_f_perp,sigma_m2_f_perp)

### mymodel3(d,k)
popt, pcov           = curve_fit(mymodel3, spacings_rw, DRW, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
k_m3_rw              = popt[0]
sigma_m3_rw          = popt[1]
popt, pcov = curve_fit(mymodel3, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m3_stat            = popt[0]
k_m3_stat_stdv       = rms_params[0]
sigma_m3_stat        = popt[1]
sigma_m3_stat_stdv   = rms_params[1]
popt, pcov           = curve_fit(mymodel3, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m3_dyn             = popt[0]
k_m3_dyn_stdv        = rms_params[0]
sigma_m3_dyn         = popt[1]
sigma_m3_dyn_stdv    = rms_params[1]
popt, pcov           = curve_fit(mymodel3, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
k_m3_f               = popt[0]
sigma_m3_f           = popt[1]
#Perp
popt, pcov              = curve_fit(mymodel3, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m3_stat_perp          = popt[0]
k_m3_stat_perp_stdv     = rms_params[0]
sigma_m3_stat_perp      = popt[1]
sigma_m3_stat_perp_stdv = rms_params[1]
popt, pcov              = curve_fit(mymodel3, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m3_dyn_perp           = popt[0]
k_m3_dyn_perp_stdv      = rms_params[0]
sigma_m3_dyn_perp       = popt[1]
sigma_m3_dyn_perp_stdv  = rms_params[1]
popt, pcov              = curve_fit(mymodel3, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=([0,-np.inf],[1,np.inf]))
k_m3_f_perp             = popt[0]
sigma_m3_f_perp         = popt[1]

# Ds
Ds_m3_rw   = mymodel3(fitds,k_m3_rw,sigma_m3_rw)
Ds_m3_stat = mymodel3(fitds,k_m3_stat,sigma_m3_stat)
Ds_m3_dyn  = mymodel3(fitds,k_m3_dyn,sigma_m3_dyn)
Ds_m3_f    = mymodel3(fitds,k_m3_f,sigma_m3_f)
Ds_m3_stat_perp = mymodel3(fitds,k_m3_stat_perp,sigma_m3_stat_perp)
Ds_m3_dyn_perp  = mymodel3(fitds,k_m3_dyn_perp,sigma_m3_dyn_perp)
Ds_m3_f_perp    = mymodel3(fitds,k_m3_f_perp,sigma_m3_f_perp)



### mymodel4(d,k,f)
popt, pcov           = curve_fit(mymodel4, spacings_rw, DRW, maxfev=thismaxfev, bounds=([0,-np.inf,-np.inf],[1,np.inf,np.inf]))
k_m4_rw              = popt[0]
f_m4_rw              = popt[1]
sigma_m4_rw          = popt[2]
popt, pcov = curve_fit(mymodel4, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=([0,-np.inf,-np.inf],[1,np.inf,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m4_stat            = popt[0]
f_m4_stat            = popt[1]
sigma_m4_stat        = popt[2]
k_m4_stat_stdv       = rms_params[0]
f_m4_stat_stdv       = rms_params[1]
sigma_m4_stat_stdv   = rms_params[2]
popt, pcov           = curve_fit(mymodel4, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=([0,-np.inf,-np.inf],[1,np.inf,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m4_dyn             = popt[0]
f_m4_dyn             = popt[1]
sigma_m4_dyn         = popt[2]
k_m4_dyn_stdv        = rms_params[0]
f_m4_dyn_stdv        = rms_params[1]
sigma_m4_dyn_stdv    = rms_params[2]
popt, pcov           = curve_fit(mymodel4, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=([0,-np.inf,-np.inf],[1,np.inf,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m4_f               = popt[0]
f_m4_f               = popt[1]
sigma_m4_f           = popt[2]
#Perp
popt, pcov              = curve_fit(mymodel4, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=([0,-np.inf,-np.inf],[1,np.inf,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m4_stat_perp          = popt[0]
f_m4_stat_perp          = popt[1]
sigma_m4_stat_perp      = popt[2]
k_m4_stat_perp_stdv     = rms_params[0]
f_m4_stat_perp_stdv     = rms_params[1]
sigma_m4_stat_perp_stdv = rms_params[2]
popt, pcov              = curve_fit(mymodel4, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=([0,-np.inf,-np.inf],[1,np.inf,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m4_dyn_perp           = popt[0]
f_m4_dyn_perp           = popt[1]
sigma_m4_dyn_perp       = popt[2]
k_m4_dyn_perp_stdv      = rms_params[0]
f_m4_dyn_perp_stdv      = rms_params[1]
sigma_m4_dyn_perp_stdv  = rms_params[2]
popt, pcov              = curve_fit(mymodel4, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=([0,-np.inf,-np.inf],[1,np.inf,np.inf]))
k_m4_f_perp             = popt[0]
f_m4_f_perp             = popt[1]
sigma_m4_f_perp         = popt[2]

# Ds
Ds_m4_rw   = mymodel4(fitds,k_m4_rw,f_m4_rw,sigma_m4_rw)
Ds_m4_stat = mymodel4(fitds,k_m4_stat,f_m4_stat,sigma_m4_stat)
Ds_m4_dyn  = mymodel4(fitds,k_m4_dyn,f_m4_dyn,sigma_m4_dyn)
Ds_m4_f    = mymodel4(fitds,k_m4_f,f_m4_f,sigma_m4_f)
Ds_m4_stat_perp = mymodel4(fitds,k_m4_stat_perp,f_m4_stat_perp,sigma_m4_stat_perp)
Ds_m4_dyn_perp  = mymodel4(fitds,k_m4_dyn_perp,f_m4_dyn_perp,sigma_m4_dyn_perp)
Ds_m4_f_perp    = mymodel4(fitds,k_m4_f_perp,f_m4_f_perp,sigma_m4_f_perp)


#'''
### mymodel5(d,k,f) # Doesn't work.
popt, pcov         = curve_fit(mymodel5, spacings_rw, DRW, maxfev=1000000, bounds=([0,-np.inf,-np.inf],[1,np.inf,np.inf]))
k_m5_rw            = popt[0]
f_m5_rw            = popt[1]
sigma_m5_rw        = popt[2]
popt, pcov         = curve_fit(mymodel5, spacings_stat, Dparallel_stat, maxfev=1000000, bounds=([0,-np.inf,-np.inf],[1,np.inf,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m5_stat          = popt[0]
f_m5_stat          = popt[1]
sigma_m5_stat      = popt[2]
k_m5_stat_stdv     = rms_params[0]
f_m5_stat_stdv     = rms_params[1]
sigma_m5_stat_stdv = rms_params[2]
popt, pcov         = curve_fit(mymodel5, spacings_dyn, Dparallel_dyn, maxfev=1000000, bounds=([0,-np.inf,-np.inf],[1,np.inf,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m5_dyn           = popt[0]
f_m5_dyn           = popt[1]
sigma_m5_dyn       = popt[2]
k_m5_dyn_stdv      = rms_params[0]
f_m5_dyn_stdv      = rms_params[1]
sigma_m5_dyn_stdv  = rms_params[2]
popt, pcov         = curve_fit(mymodel5, spacings_f, Dparallel_f, maxfev=1000000, bounds=([0,-np.inf,-np.inf],[1,np.inf,np.inf]))
k_m5_f             = popt[0]
f_m5_f             = popt[1]
sigma_m5_f         = popt[2]
#Perp
popt, pcov              = curve_fit(mymodel5, spacings_stat, Dzs_stat, maxfev=1000000, bounds=([0,-np.inf,-np.inf],[1,np.inf,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m5_stat_perp          = popt[0]
f_m5_stat_perp          = popt[1]
sigma_m5_stat_perp      = popt[2]
k_m5_stat_perp_stdv     = rms_params[0]
f_m5_stat_perp_stdv     = rms_params[1]
sigma_m5_stat_perp_stdv = rms_params[2]
popt, pcov              = curve_fit(mymodel5, spacings_dyn, Dzs_dyn, maxfev=1000000, bounds=([0,-np.inf,-np.inf],[1,np.inf,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m5_dyn_perp           = popt[0]
f_m5_dyn_perp           = popt[1]
sigma_m5_dyn_perp       = popt[2]
k_m5_dyn_perp_stdv      = rms_params[0]
f_m5_dyn_perp_stdv      = rms_params[1]
sigma_m5_dyn_perp_stdv  = rms_params[2]
popt, pcov              = curve_fit(mymodel5, spacings_f, Dzs_f, maxfev=1000000, bounds=([0,-np.inf,-np.inf],[1,np.inf,np.inf]))
k_m5_f_perp             = popt[0]
f_m5_f_perp             = popt[1]
sigma_m5_f_perp         = popt[2]

# Ds
Ds_m5_rw   = mymodel5(fitds,k_m5_rw,f_m5_rw,sigma_m5_rw)
Ds_m5_stat = mymodel5(fitds,k_m5_stat,f_m5_stat,sigma_m5_stat)
Ds_m5_dyn  = mymodel5(fitds,k_m5_dyn,f_m5_dyn,sigma_m5_dyn)
Ds_m5_f    = mymodel5(fitds,k_m5_f,f_m5_f,sigma_m5_f)
Ds_m5_stat_perp = mymodel5(fitds,k_m5_stat_perp,f_m5_stat_perp,sigma_m5_stat_perp)
Ds_m5_dyn_perp  = mymodel5(fitds,k_m5_dyn_perp,f_m5_dyn_perp,sigma_m5_dyn_perp)
Ds_m5_f_perp    = mymodel5(fitds,k_m5_f_perp,f_m5_f_perp,sigma_m5_f_perp)
#'''

### mymodel6(d,k)
popt, pcov         = curve_fit(mymodel6, spacings_rw, DRW, maxfev=thismaxfev)#, bounds=(0,np.inf))
k_m6_rw            = popt[0]
sigma_m6_rw        = popt[1]
popt, pcov         = curve_fit(mymodel6, spacings_stat, Dparallel_stat, maxfev=thismaxfev)#, bounds=(0,np.inf))
rms_params = np.sqrt(np.diag(pcov))
k_m6_stat          = popt[0]
k_m6_stat_stdv     = rms_params[0]
sigma_m6_stat      = popt[1]
sigma_m6_stat_stdv = rms_params[1]
popt, pcov         = curve_fit(mymodel6, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev)#, bounds=(0,np.inf))
rms_params = np.sqrt(np.diag(pcov))
k_m6_dyn           = popt[0]
k_m6_dyn_stdv      = rms_params[0]
sigma_m6_dyn       = popt[1]
sigma_m6_dyn_stdv  = rms_params[1]
popt, pcov         = curve_fit(mymodel6, spacings_f, Dparallel_f, maxfev=thismaxfev)#, bounds=(0,np.inf))
k_m6_f             = popt[0]
sigma_m6_f         = popt[1]
#Perp
popt, pcov              = curve_fit(mymodel6, spacings_stat, Dzs_stat, maxfev=thismaxfev)#, bounds=(0,np.inf))
rms_params = np.sqrt(np.diag(pcov))
k_m6_stat_perp          = popt[0]
k_m6_stat_perp_stdv     = rms_params[0]
sigma_m6_stat_perp      = popt[1]
sigma_m6_stat_perp_stdv = rms_params[1]
popt, pcov              = curve_fit(mymodel6, spacings_dyn, Dzs_dyn, maxfev=thismaxfev)#, bounds=(0,np.inf))
rms_params = np.sqrt(np.diag(pcov))
k_m6_dyn_perp           = popt[0]
k_m6_dyn_perp_stdv      = rms_params[0]
sigma_m6_dyn_perp       = popt[1]
sigma_m6_dyn_perp_stdv  = rms_params[1]
#print('Stdv:', np.sqrt(np.diag(pcov)))
popt, pcov              = curve_fit(mymodel6, spacings_f, Dzs_f, maxfev=thismaxfev)#, bounds=(0,np.inf))
k_m6_f_perp             = popt[0]
sigma_m6_f_perp         = popt[1]

# Ds
Ds_m6_rw   = mymodel6(fitds,k_m6_rw,sigma_m6_rw)
Ds_m6_stat = mymodel6(fitds,k_m6_stat,sigma_m6_stat)
Ds_m6_dyn  = mymodel6(fitds,k_m6_dyn,sigma_m6_dyn)
Ds_m6_f    = mymodel6(fitds,k_m6_f,sigma_m6_f)
Ds_m6_stat_perp = mymodel6(fitds,k_m6_stat_perp,sigma_m6_stat_perp)
Ds_m6_dyn_perp  = mymodel6(fitds,k_m6_dyn_perp,sigma_m6_dyn_perp)
Ds_m6_f_perp    = mymodel6(fitds,k_m6_f_perp,sigma_m6_f_perp)

### mymodel7(d,k)
popt, pcov         = curve_fit(mymodel7, spacings_rw, DRW, maxfev=thismaxfev, bounds=([0,0],[1,np.inf]))
k_m7_rw            = popt[0]
sigma_m7_rw        = popt[1]
popt, pcov         = curve_fit(mymodel7, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=([0,0],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m7_stat          = popt[0]
k_m7_stat_stdv     = rms_params[0]
sigma_m7_stat      = popt[1]
sigma_m7_stat_stdv = rms_params[1]
popt, pcov         = curve_fit(mymodel7, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=([0,0],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m7_dyn           = popt[0]
k_m7_dyn_stdv      = rms_params[0]
sigma_m7_dyn       = popt[1]
sigma_m7_dyn_stdv  = rms_params[1]
popt, pcov         = curve_fit(mymodel7, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=([0,0],[1,np.inf]))
k_m7_f             = popt[0]
sigma_m7_f         = popt[1]
#Perp
popt, pcov              = curve_fit(mymodel7, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=([0,0],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m7_stat_perp          = popt[0]
k_m7_stat_perp_stdv     = rms_params[0]
sigma_m7_stat_perp      = popt[1]
sigma_m7_stat_perp_stdv = rms_params[1]
popt, pcov              = curve_fit(mymodel7, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=([0,0],[1,np.inf]))
rms_params = np.sqrt(np.diag(pcov))
k_m7_dyn_perp           = popt[0]
k_m7_dyn_perp_stdv      = rms_params[0]
sigma_m7_dyn_perp       = popt[1]
sigma_m7_dyn_perp_stdv  = rms_params[1]
popt, pcov              = curve_fit(mymodel7, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=([0,0],[1,np.inf]))
k_m7_f_perp             = popt[0]
sigma_m7_f_perp         = popt[1]

# Ds
Ds_m7_rw   = mymodel7(fitds,k_m7_rw,sigma_m7_rw)
Ds_m7_stat = mymodel7(fitds,k_m7_stat,sigma_m7_stat)
Ds_m7_dyn  = mymodel7(fitds,k_m7_dyn,sigma_m7_dyn)
Ds_m7_f    = mymodel7(fitds,k_m7_f,sigma_m7_f)
Ds_m7_stat_perp = mymodel7(fitds,k_m7_stat_perp,sigma_m7_stat_perp)
Ds_m7_dyn_perp  = mymodel7(fitds,k_m7_dyn_perp,sigma_m7_dyn_perp)
Ds_m7_f_perp    = mymodel7(fitds,k_m7_f_perp,sigma_m7_f_perp)


### powerlaw1(d,k,f)
popt, pcov       = curve_fit(powerlaw1, spacings_rw, DRW, maxfev=thismaxfev)#, bounds=(0,1))
k_p1_rw          = popt[0]
l_p1_rw          = popt[1]
n_p1_rw          = popt[2]
popt, pcov       = curve_fit(powerlaw1, spacings_stat, Dparallel_stat, maxfev=thismaxfev)#, bounds=(0,1))
rms_params = np.sqrt(np.diag(pcov))
k_p1_stat        = popt[0]
l_p1_stat        = popt[1]
n_p1_stat        = popt[2]
k_p1_stat_stdv   = rms_params[0]
l_p1_stat_stdv   = rms_params[1]
n_p1_stat_stdv   = rms_params[2]
popt, pcov       = curve_fit(powerlaw1, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev)#, bounds=(0,1))
rms_params = np.sqrt(np.diag(pcov))
k_p1_dyn         = popt[0]
l_p1_dyn         = popt[1]
n_p1_dyn         = popt[2]
k_p1_dyn_stdv    = rms_params[0]
l_p1_dyn_stdv    = rms_params[1]
n_p1_dyn_stdv    = rms_params[2]
popt, pcov       = curve_fit(powerlaw1, spacings_f, Dparallel_f, maxfev=thismaxfev)#, bounds=(0,1))
k_p1_f           = popt[0]
l_p1_f           = popt[1]
n_p1_f           = popt[2]
#Perp
popt, pcov          = curve_fit(powerlaw1, spacings_stat, Dzs_stat, maxfev=thismaxfev)#, bounds=(0,1))
rms_params = np.sqrt(np.diag(pcov))
k_p1_stat_perp      = popt[0]
l_p1_stat_perp      = popt[1]
n_p1_stat_perp      = popt[2]
k_p1_stat_perp_stdv = rms_params[0]
l_p1_stat_perp_stdv = rms_params[1]
n_p1_stat_perp_stdv = rms_params[2]
popt, pcov          = curve_fit(powerlaw1, spacings_dyn, Dzs_dyn, maxfev=thismaxfev)#, bounds=(0,1))
rms_params = np.sqrt(np.diag(pcov))
k_p1_dyn_perp       = popt[0]
l_p1_dyn_perp       = popt[1]
n_p1_dyn_perp       = popt[2]
k_p1_dyn_perp_stdv  = rms_params[0]
l_p1_dyn_perp_stdv  = rms_params[1]
n_p1_dyn_perp_stdv  = rms_params[2]
popt, pcov          = curve_fit(powerlaw1, spacings_f, Dzs_f, maxfev=thismaxfev)#, bounds=(0,1))
k_p1_f_perp         = popt[0]
l_p1_f_perp         = popt[1]
n_p1_f_perp         = popt[2]

# Ds
Ds_p1_rw   = powerlaw1(fitds,k_p1_rw,l_p1_rw,n_p1_rw)
Ds_p1_stat = powerlaw1(fitds,k_p1_stat,l_p1_stat,n_p1_stat)
Ds_p1_dyn  = powerlaw1(fitds,k_p1_dyn,l_p1_dyn,n_p1_dyn)
Ds_p1_f    = powerlaw1(fitds,k_p1_f,l_p1_f,n_p1_f)
Ds_p1_stat_perp = powerlaw1(fitds,k_p1_stat_perp,l_p1_stat_perp,n_p1_stat_perp)
Ds_p1_dyn_perp  = powerlaw1(fitds,k_p1_dyn_perp,l_p1_dyn_perp,n_p1_dyn_perp)
Ds_p1_f_perp    = powerlaw1(fitds,k_p1_f_perp,l_p1_f_perp,n_p1_f_perp)


#### Plotting:
##params = {'mathtext.default': 'regular' }  # Ditch this.   
###plt.rcParams.update(params)

plt.figure(figsize=(6,5))#(14,5))
ax = plt.subplot(111)
ax.plot(spacings_dyn, Dparallel_dyn, color='g', label=r'$D_\parallel$, dyn.')
ax.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
#ax.plot(dg,DRW, color='b', label='Random walk, fixed geom., no refl.')
#ax.plot(fitds,Ds_op_dyn, '--', label=r'Ordered packings')
#ax.plot(fitds,Ds_hr_dyn, '-.', label=r'Hyperbola of revolution') # 2
#ax.plot(fitds,Ds_nm_dyn, '--', label=r'Not monosized spheres')  ####
#ax.plot(fitds,Ds_ps_dyn, '--', label=r'Monodisperse sphere packing') ####
#ax.plot(fitds,Ds_os_dyn, '--', label=r'Overlapping spheres') # 5 ####
#ax.plot(fitds,Ds_oc_dyn, '--', label=r'Overlapping cylinders')  ####
ax.plot(fitds,Ds_hc_dyn, ':', label=r'Heterogeneous catalyst') 
#ax.plot(fitds,Ds_m1_dyn, '-.', label=r'Custom model 1')
#ax.plot(fitds,Ds_m2_dyn, '-.', label=r'Custom model 2')
#ax.plot(fitds,Ds_m3_dyn, '-.', label=r'Custom model 3')
line1, = ax.plot(fitds,Ds_rm_dyn, '--,', label=r'Cation-exchange resin membrane') 
line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
ax.plot(fitds,Ds_m4_dyn, '-.', color='rebeccapurple', label=r'Custom model')#'$\tau^2=[1+v_p/(1-kv_p^f)]^2$')#
#ax.plot(fitds,Ds_m5_dyn, '-.', label=r'Custom model 5')
#ax.plot(fitds,Ds_m6_dyn, '-.', label=r'Custom model 6')
#ax.plot(fitds,Ds_m7_dyn, '-.', label=r'Custom model 7')
#ax.plot(fitds,Ds_p1_dyn, '-.', label=r'Power law 1')
# Add more plotting! # Need to beware of negative porosity for some!
plt.xlabel(r'$d$ (nm)', fontsize=12)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=12)
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(plotname_d_dyn)


## Plotting:# Need fancy plotting. Find code and redo fancy plotting.
plt.figure(figsize=(6,5))#(14,5))
ax = plt.subplot(111)
ax.plot(spacings_stat, Dparallel_stat, color='limegreen',label=r'$D_\parallel$, stat.')
ax.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen', alpha=0.2)
#ax.plot(fitds,Ds_op_stat, '--', label=r'Ordered packings') # I have most values for the dynamic brush
#ax.plot(fitds,Ds_hr_stat, '-.', label=r'Hyperbola of revolution') # 2
#ax.plot(fitds,Ds_nm_stat, '--', label=r'Not monosized spheres')  ####
#ax.plot(fitds,Ds_ps_stat, '--', label=r'Monodisperse sphere packing') ####
#ax.plot(fitds,Ds_os_stat, '--', label=r'Overlapping spheres') # 5 ####
#ax.plot(fitds,Ds_oc_stat, '--', label=r'Overlapping cylinders')  ####
ax.plot(fitds,Ds_hc_stat, ':', label=r'Heterogeneous catalyst') 
#ax.plot(fitds,Ds_m1_stat, '-.', label=r'Custom model 1')
#ax.plot(fitds,Ds_m2_stat, '-.', label=r'Custom model 2')
#ax.plot(fitds,Ds_m3_stat, '-.', label=r'Custom model 3')
line1, = ax.plot(fitds,Ds_rm_stat, '--,', label=r'Cation-exchange resin membrane') 
line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
ax.plot(fitds,Ds_m4_stat, '-.', color='rebeccapurple', label=r'Custom model')#'$\tau^2=[1+v_p/(1-kv_p^f)]^2$')#
#ax.plot(fitds,Ds_m5_stat, '-.', label=r'Custom model 5')
#ax.plot(fitds,Ds_m6_stat, '-.', label=r'Custom model 6')
#ax.plot(fitds,Ds_m7_stat, '-.', label=r'Custom model 7')
#ax.plot(fitds,Ds_p1_stat, '-.', label=r'Power law 1')
# Add more plotting! # Need to beware of negative porosity for some!
plt.xlabel(r'$d$ (nm)')
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=12)
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(plotname_d_stat)




## Plotting:# Need fancy plotting. Find code and redo fancy plotting.
plt.figure(figsize=(6,5))#(14,5))
ax = plt.subplot(111)
ax.plot(spacings_f, Dparallel_f, color='springgreen', label=r'$D_\parallel$, straight')
ax.fill_between(spacings_f, Dparallel_f+Dparallel_stdv_f, Dparallel_f-Dparallel_stdv_f, facecolor='springgreen', alpha=0.2)
#ax.plot(fitds,Ds_op_f, '--', label=r'Ordered packings') # I have most values for the dynamic brush
#ax.plot(fitds,Ds_hr_f, '-.', label=r'Hyperbola of revolution') # 2
#ax.plot(fitds,Ds_nm_f, '--', label=r'Not monosized spheres')  ####
#ax.plot(fitds,Ds_ps_f, '--', label=r'Monodisperse sphere packing') ####
#ax.plot(fitds,Ds_os_f, '--', label=r'Overlapping spheres') # 5 ####
#ax.plot(fitds,Ds_oc_f, '--', label=r'Overlapping cylinders')  ####
ax.plot(fitds,Ds_hc_f, ':', label=r'Heterogeneous catalyst')
#ax.plot(fitds,Ds_m1_f, '-.', label=r'Custom model 1')
#ax.plot(fitds,Ds_m2_f, '-.', label=r'Custom model 2')
#ax.plot(fitds,Ds_m3_f, '-.', label=r'Custom model 3')
line1, = ax.plot(fitds,Ds_rm_f, '--,', label=r'Cation-exchange resin membrane') 
line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
ax.plot(fitds,Ds_m4_f, '-.', color='rebeccapurple', label=r'Custom model')#'$\tau^2=[1+v_p/(1-kv_p^f)]^2$')#
#ax.plot(fitds,Ds_m5_f, '-.', label=r'Custom model 5')
#ax.plot(fitds,Ds_m6_f, '-.', label=r'Custom model 6')
#ax.plot(fitds,Ds_m7_f, '-.', label=r'Custom model 7')
#ax.plot(fitds,Ds_p1_f, '-.', label=r'Power law 1')
# Add more plotting! # Need to beware of negative porosity for some!
plt.xlabel(r'$d$ (nm)', fontsize=12)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=12)
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(plotname_d_f)

#... Do the same, men for perp.
plt.figure(figsize=(6,5))#(14,5))
ax = plt.subplot(111)
ax.plot(spacings_dyn, Dzs_dyn, color='g', label=r'$D_\perp$, dyn.')
ax.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='g', alpha=0.2)
#ax.plot(fitds,Ds_op_dyn_perp, '--', label=r'Ordered packings')
#ax.plot(fitds,Ds_hr_dyn_perp, '-.', label=r'Hyperbola of revolution') # 2
#ax.plot(fitds,Ds_nm_dyn_perp, '--', label=r'Not monosized spheres')  ####
#ax.plot(fitds,Ds_ps_dyn_perp, '--', label=r'Monodisperse sphere packing') ####
#ax.plot(fitds,Ds_os_dyn_perp, '--', label=r'Overlapping spheres') # 5 ####
#ax.plot(fitds,Ds_oc_dyn_perp, '--', label=r'Overlapping cylinders')  ####
ax.plot(fitds,Ds_hc_dyn_perp, ':', label=r'Heterogeneous catalyst') 
#ax.plot(fitds,Ds_m1_dyn_perp, '-.', label=r'Custom model 1')
#ax.plot(fitds,Ds_m2_dyn_perp, '-.', label=r'Custom model 2')
#ax.plot(fitds,Ds_m3_dyn_perp, '-.', label=r'Custom model 3')
line1, = ax.plot(fitds,Ds_rm_dyn_perp, '--,', label=r'Cation-exchange resin membrane') 
line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
ax.plot(fitds,Ds_m4_dyn_perp, '-.', color='rebeccapurple', label=r'Custom model')#'$\tau^2=[1+v_p/(1-kv_p^f)]^2$')#
#ax.plot(fitds,Ds_m5_dyn_perp, '-.', label=r'Custom model 5')
#ax.plot(fitds,Ds_m6_dyn_perp, '-.', label=r'Custom model 6')
#ax.plot(fitds,Ds_m7_dyn_perp, '-.', label=r'Custom model 7')
#ax.plot(fitds,Ds_p1_dyn_perp, '-.', label=r'Power law 1')
# Add more plotting! # Need to beware of negative porosity for some!
plt.xlabel(r'$d$ (nm)', fontsize=12)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=12)
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(plotname_d_dyn_perp)


plt.figure(figsize=(6,5))#(14,5))
ax = plt.subplot(111)
ax.plot(spacings_stat, Dzs_stat, color='limegreen', label=r'$D_\perp$, stat.')
ax.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='limegreen', alpha=0.2)
#ax.plot(fitds,Ds_op_stat_perp, '--', label=r'Ordered packings') # I have most values for the dynamic brush
#ax.plot(fitds,Ds_hr_stat_perp, '-.', label=r'Hyperbola of revolution') # 2
#ax.plot(fitds,Ds_nm_stat_perp, '--', label=r'Not monosized spheres')  ####
#ax.plot(fitds,Ds_ps_stat_perp, '--', label=r'Monodisperse sphere packing') ####
#ax.plot(fitds,Ds_os_stat_perp, '--', label=r'Overlapping spheres') # 5 ####
#ax.plot(fitds,Ds_oc_stat_perp, '--', label=r'Overlapping cylinders')  ####
ax.plot(fitds,Ds_hc_stat_perp, ':', label=r'Heterogeneous catalyst') 
#ax.plot(fitds,Ds_m1_stat_perp, '-.', label=r'Custom model 1')
#ax.plot(fitds,Ds_m2_stat_perp, '-.', label=r'Custom model 2')
#ax.plot(fitds,Ds_m3_stat_perp, '-.', label=r'Custom model 3')
line1, = ax.plot(fitds,Ds_rm_stat_perp, '--,', label=r'Cation-exchange resin membrane') 
line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
ax.plot(fitds,Ds_m4_stat_perp, '-.', color='rebeccapurple', label=r'Custom model')#'$\tau^2=[1+v_p/(1-kv_p^f)]^2$')#
#ax.plot(fitds,Ds_m5_stat_perp, '-.', label=r'Custom model 5')
#ax.plot(fitds,Ds_m6_stat_perp, '-.', label=r'Custom model 6')
#ax.plot(fitds,Ds_m7_stat_perp, '-.', label=r'Custom model 7')
#ax.plot(fitds,Ds_p1_stat_perp, '-.', label=r'Power law 1')
# Add more plotting! # Need to beware of negative porosity for some!
plt.xlabel(r'$d$ (nm)', fontsize=12)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=12)
#plt.title(r'$D/D_{\mathregular{bulk}}$ vs $d$')
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(plotname_d_stat_perp)

plt.figure(figsize=(6,5))#(14,5))
ax = plt.subplot(111)
ax.errorbar(spacings_f, Dzs_f, color='springgreen', label=r'$D_\perp$, straight')
ax.fill_between(spacings_f, Dzs_f+Dzs_stdv_f, Dzs_f-Dzs_stdv_f, facecolor='g', alpha=0.2)
#ax.plot(fitds,Ds_op_f_perp, '--', label=r'Ordered packings') # I have most values for the dynamic brush
#ax.plot(fitds,Ds_hr_f_perp, '-.', label=r'Hyperbola of revolution') # 2
#ax.plot(fitds,Ds_nm_f_perp, '--', label=r'Not monosized spheres')  ####
#ax.plot(fitds,Ds_ps_f_perp, '--', label=r'Monodisperse sphere packing') ####
#ax.plot(fitds,Ds_os_f_perp, '--', label=r'Overlapping spheres') # 5 ####
#ax.plot(fitds,Ds_oc_f_perp, '--', label=r'Overlapping cylinders')  ####
ax.plot(fitds,Ds_hc_f_perp, ':', label=r'Heterogeneous catalyst')
#ax.plot(fitds,Ds_m1_f_perp, '-.', label=r'Custom model 1')
#ax.plot(fitds,Ds_m2_f_perp, '-.', label=r'Custom model 2')
#ax.plot(fitds,Ds_m3_f_perp, '-.', label=r'Custom model 3')
line1, = ax.plot(fitds,Ds_rm_f_perp, '--,', label=r'Cation-exchange resin membrane') 
line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
ax.plot(fitds,Ds_m4_f_perp, '-.', color='rebeccapurple', label=r'Custom model')#'$\tau^2=[1+v_p/(1-kv_p^f)]^2$')#
#ax.plot(fitds,Ds_m5_f_perp, '-.', label=r'Custom model 5')
#ax.plot(fitds,Ds_m6_f_perp, '-.', label=r'Custom model 6')
#ax.plot(fitds,Ds_m7_f_perp, '-.', label=r'Custom model 7')
#ax.plot(fitds,Ds_p1_f_perp, '-.', label=r'Power law 1')
# Add more plotting! # Need to beware of negative porosity for some!
plt.xlabel(r'$d$ (nm)', fontsize=12)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=12)
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(plotname_d_f_perp)

print('\n\nModel 1')
print('k stat, par:', k_m1_stat, '+/-', k_m1_stat_stdv)
print('k stat, perp:', k_m1_stat_perp, '+/-', k_m1_stat_perp_stdv)
print('k dyn, par:', k_m1_dyn, '+/-', k_m1_dyn_stdv)
print('k dyn, perp:', k_m1_dyn_perp, '+/-', k_m1_dyn_perp_stdv)

print('sigma stat, par:', sigma_m1_stat, '+/-', sigma_m1_stat_stdv)
print('sigma stat, perp:', sigma_m1_stat_perp, '+/-', sigma_m1_stat_perp_stdv)
print('sigma dyn, par:', sigma_m1_dyn, '+/-', sigma_m1_dyn_stdv)
print('sigma dyn, perp:', sigma_m1_dyn_perp, '+/-', sigma_m1_dyn_perp_stdv)

print('\n\nModel 2')
print('k stat, par:', k_m2_stat, '+/-', k_m2_stat_stdv)
print('k stat, perp:', k_m2_stat_perp, '+/-', k_m2_stat_perp_stdv)
print('k dyn, par:', k_m2_dyn, '+/-', k_m2_dyn_stdv)
print('k dyn, perp:', k_m2_dyn_perp, '+/-', k_m2_dyn_perp_stdv)

print('sigma stat, par:', sigma_m2_stat, '+/-', sigma_m2_stat_stdv)
print('sigma stat, perp:', sigma_m2_stat_perp, '+/-', sigma_m2_stat_perp_stdv)
print('sigma dyn, par:', sigma_m2_dyn, '+/-', sigma_m2_dyn_stdv)
print('sigma dyn, perp:', sigma_m2_dyn_perp, '+/-', sigma_m2_dyn_perp_stdv)


# Model 3
print('\n\nModel 3')
print('k stat, par:', k_m3_stat, '+/-', k_m3_stat_stdv)
print('k stat, perp:', k_m3_stat_perp, '+/-', k_m3_stat_perp_stdv)
print('k dyn, par:', k_m3_dyn, '+/-', k_m3_dyn_stdv)
print('k dyn, perp:', k_m3_dyn_perp, '+/-', k_m3_dyn_perp_stdv)


print('sigma stat, par:', sigma_m3_stat, '+/-', sigma_m3_stat_stdv)
print('sigma stat, perp:', sigma_m3_stat_perp, '+/-', sigma_m3_stat_perp_stdv)
print('sigma dyn, par:', sigma_m3_dyn, '+/-', sigma_m3_dyn_stdv)
print('sigma dyn, perp:', sigma_m3_dyn_perp, '+/-', sigma_m3_dyn_perp_stdv)


# Model 4
print('\n\nModel 4')
print('k stat, par:', k_m4_stat, '+/-', k_m4_stat_stdv,'; f stat, par:', f_m4_stat, '+/-', f_m4_stat_stdv)
print('k stat, perp:', k_m4_stat_perp, '+/-', k_m4_stat_perp_stdv,'; f stat, perp:', f_m4_stat_perp, '+/-', f_m4_stat_perp_stdv)
print('k dyn, par:', k_m4_dyn, '+/-', k_m4_dyn_stdv,'; f dyn, par:', f_m4_dyn, '+/-', f_m4_dyn_stdv)
print('k dyn, perp:', k_m4_dyn_perp, '+/-', k_m4_dyn_perp_stdv,'; f dyn, perp:', f_m4_dyn_perp, '+/-', f_m4_dyn_perp_stdv)

print('sigma stat, par:', sigma_m4_stat, '+/-', sigma_m4_stat_stdv)
print('sigma stat, perp:', sigma_m4_stat_perp, '+/-', sigma_m4_stat_perp_stdv)
print('sigma dyn, par:', sigma_m4_dyn, '+/-', sigma_m4_dyn_stdv)
print('sigma dyn, perp:', sigma_m4_dyn_perp, '+/-', sigma_m4_dyn_perp_stdv)

#'''
# Model 5
print('\n\nModel 5')
print('k stat, par:', k_m5_stat, '+/-', k_m5_stat_stdv,'; f stat, par:', f_m5_stat, '+/-', f_m5_stat_stdv)
print('k stat, perp:', k_m5_stat_perp, '+/-', k_m5_stat_perp_stdv,'; f stat, perp:', f_m5_stat_perp, '+/-', f_m5_stat_perp_stdv)
print('k dyn, par:', k_m5_dyn, '+/-', k_m5_dyn_stdv,'; f dyn, par:', f_m5_dyn, '+/-', f_m5_dyn_stdv)
print('k dyn, perp:', k_m5_dyn_perp, '+/-', k_m5_dyn_perp_stdv,'; f dyn, perp:', f_m5_dyn_perp, '+/-', f_m5_dyn_perp_stdv)

print('sigma stat, par:', sigma_m5_stat, '+/-', sigma_m5_stat_stdv)
print('sigma stat, perp:', sigma_m5_stat_perp, '+/-', sigma_m5_stat_perp_stdv)
print('sigma dyn, par:', sigma_m5_dyn, '+/-', sigma_m5_dyn_stdv)
print('sigma dyn, perp:', sigma_m5_dyn_perp, '+/-', sigma_m5_dyn_perp_stdv)

#'''

# Model 6
print('\n\nModel 6')
print('k stat, par:', k_m6_stat, '+/-', k_m6_stat_stdv)
print('k stat, perp:', k_m6_stat_perp, '+/-', k_m6_stat_perp_stdv)
print('k dyn, par:', k_m6_dyn, '+/-', k_m6_dyn_stdv)
print('k dyn, perp:', k_m6_dyn_perp, '+/-', k_m6_dyn_perp_stdv)

print('sigma stat, par:', sigma_m6_stat, '+/-', sigma_m6_stat_stdv)
print('sigma stat, perp:', sigma_m6_stat_perp, '+/-', sigma_m6_stat_perp_stdv)
print('sigma dyn, par:', sigma_m6_dyn, '+/-', sigma_m6_dyn_stdv)
print('sigma dyn, perp:', sigma_m6_dyn_perp, '+/-', sigma_m6_dyn_perp_stdv)


# Model 7
print('\n\nModel 7')
print('k stat, par:', k_m7_stat, '+/-', k_m7_stat_stdv)
print('k stat, perp:', k_m7_stat_perp, '+/-', k_m7_stat_perp_stdv)
print('k dyn, par:', k_m7_dyn, '+/-', k_m7_dyn_stdv)
print('k dyn, perp:', k_m7_dyn_perp, '+/-', k_m7_dyn_perp_stdv)

print('sigma stat, par:', sigma_m7_stat, '+/-', sigma_m7_stat_stdv)
print('sigma stat, perp:', sigma_m7_stat_perp, '+/-', sigma_m7_stat_perp_stdv)
print('sigma dyn, par:', sigma_m7_dyn, '+/-', sigma_m7_dyn_stdv)
print('sigma dyn, perp:', sigma_m7_dyn_perp, '+/-', sigma_m7_dyn_perp_stdv)
print('  ')

# Shen & Chen-models:

print('Ordered packings:')
print('Sigma, stat, par:', sigma_op_stat, '+/-', sigma_op_stat_rms)
print('Sigma, stat, perp:', sigma_op_stat_perp, '+/-', sigma_op_stat_perp_rms)
print('Sigma, dyn, par:', sigma_op_dyn, '+/-', sigma_op_dyn_rms)
print('Sigma, dyn, perp:', sigma_op_dyn_perp, '+/-', sigma_op_dyn_perp_rms)
print('   ')

print('Hyperbola of revolution:')
print('Sigma, stat, par:', sigma_hr_stat, '+/-', sigma_hr_stat_rms)
print('Sigma, stat, perp:', sigma_hr_stat_perp, '+/-', sigma_hr_stat_perp_rms)
print('Sigma, dyn, par:', sigma_hr_dyn, '+/-', sigma_hr_dyn_rms)
print('Sigma, dyn, perp:', sigma_hr_dyn_perp, '+/-', sigma_hr_dyn_perp_rms)
print('   ')

print('Heterogeneous catalyst:')
print('Sigma, stat, par:', sigma_hc_stat, '+/-', sigma_hc_stat_rms)
print('Sigma, stat, perp:', sigma_hc_stat_perp, '+/-', sigma_hc_stat_perp_rms)
print('Sigma, dyn, par:', sigma_hc_dyn, '+/-', sigma_hc_dyn_rms)
print('Sigma, dyn, perp:', sigma_hc_dyn_perp, '+/-', sigma_hc_dyn_perp_rms)
print('   ')

print('Cation-exchange resin membrane:')
print('Sigma, stat, par:', sigma_rm_stat, '+/-', sigma_rm_stat_rms)
print('Sigma, stat, perp:', sigma_rm_stat_perp, '+/-', sigma_rm_stat_perp_rms)
print('Sigma, dyn, par:', sigma_rm_dyn, '+/-', sigma_rm_dyn_rms)
print('Sigma, dyn, perp:', sigma_rm_dyn_perp, '+/-', sigma_rm_dyn_perp_rms)
print('   ')


print('-----------------------')
print('Ordered packings:')
#
phi = giveporosity(fitds,sigma_op_dyn)
minphi = min(phi); maxphi = max(phi)
print('op, dyn, par. min(phi):',minphi,', max(phi):',maxphi)
#
phi = giveporosity(fitds,sigma_op_stat)
minphi = min(phi); maxphi = max(phi)
print('op, stat, par. min(phi):',minphi,', max(phi):',maxphi)
#
phi = giveporosity(fitds,sigma_op_dyn_perp)
minphi = min(phi); maxphi = max(phi)
print('op, dyn, perp. min(phi):',minphi,', max(phi):',maxphi)
#
phi = giveporosity(fitds,sigma_op_stat_perp)
minphi = min(phi); maxphi = max(phi)
print('op, stat, perp. min(phi):',minphi,', max(phi):',maxphi)
print('-----------------------')
print('Hyperbola of revolution:')
phi = giveporosity(fitds,sigma_hr_dyn)
minphi = min(phi); maxphi = max(phi)
print('hr, dyn, par. min(phi):',minphi,', max(phi):',maxphi)
#
phi = giveporosity(fitds,sigma_hr_stat)
minphi = min(phi); maxphi = max(phi)
print('hr, stat, par. min(phi):',minphi,', max(phi):',maxphi)
#
phi = giveporosity(fitds,sigma_hr_dyn_perp)
minphi = min(phi); maxphi = max(phi)
print('hr, dyn, perp. min(phi):',minphi,', max(phi):',maxphi)
#
phi = giveporosity(fitds,sigma_hr_stat_perp)
minphi = min(phi); maxphi = max(phi)
print('hr, stat, perp. min(phi):',minphi,', max(phi):',maxphi)
#
print('-----------------------')
print('Heterogeneous catalyst:')
phi = giveporosity(fitds,sigma_hc_dyn)
minphi = min(phi); maxphi = max(phi)
print('hc, dyn, par. min(phi):',minphi,', max(phi):',maxphi)
#
phi = giveporosity(fitds,sigma_hc_stat)
minphi = min(phi); maxphi = max(phi)
print('hc, stat, par. min(phi):',minphi,', max(phi):',maxphi)
#
phi = giveporosity(fitds,sigma_hc_dyn_perp)
minphi = min(phi); maxphi = max(phi)
print('hc, dyn, perp. min(phi):',minphi,', max(phi):',maxphi)
#
phi = giveporosity(fitds,sigma_hc_stat_perp)
minphi = min(phi); maxphi = max(phi)
print('hc, stat, perp. min(phi):',minphi,', max(phi):',maxphi)
print('-----------------------')
print('Cation-exchange resin membrane:')
phi = giveporosity(fitds,sigma_rm_dyn)
minphi = min(phi); maxphi = max(phi)
print('rm, dyn, par. min(phi):',minphi,', max(phi):',maxphi)
#
phi = giveporosity(fitds,sigma_rm_stat)
minphi = min(phi); maxphi = max(phi)
print('rm, stat, par. min(phi):',minphi,', max(phi):',maxphi)
#
phi = giveporosity(fitds,sigma_rm_dyn_perp)
minphi = min(phi); maxphi = max(phi)
print('rm, dyn, perp. min(phi):',minphi,', max(phi):',maxphi)
#
phi = giveporosity(fitds,sigma_rm_stat_perp)
minphi = min(phi); maxphi = max(phi)
print('rm, stat, perp. min(phi):',minphi,', max(phi):',maxphi)
print('-----------------------')


print('-----------------------')
print('Ordered packings:')
#
phi = giveporosity(fitds,sigma_op_dyn)
passd = find_d_for_first_positive(phi,fitds)
print('op, dyn, par. passd:',passd)
#
phi = giveporosity(fitds,sigma_op_stat)
passd = find_d_for_first_positive(phi,fitds)
print('op, stat, par. passd:',passd)
#
phi = giveporosity(fitds,sigma_op_dyn_perp)
passd = find_d_for_first_positive(phi,fitds)
print('op, dyn, perp. passd:',passd)
#
phi = giveporosity(fitds,sigma_op_stat_perp)
passd = find_d_for_first_positive(phi,fitds)
print('op, dyn, perp. passd:',passd)
print('-----------------------')
print('Hyperbola of revolution:')
phi = giveporosity(fitds,sigma_hr_dyn)
passd = find_d_for_first_positive(phi,fitds)
print('hr, dyn, par. passd:',passd)
#
phi = giveporosity(fitds,sigma_hr_stat)
passd = find_d_for_first_positive(phi,fitds)
print('hr, stat, par. passd:',passd)
#
phi = giveporosity(fitds,sigma_hr_dyn_perp)
passd = find_d_for_first_positive(phi,fitds)
print('hr, dyn, perp. passd:',passd)
#
phi = giveporosity(fitds,sigma_hr_stat_perp)
passd = find_d_for_first_positive(phi,fitds)
print('hr, stat, perp. passd:',passd)
#
print('-----------------------')
print('Heterogeneous catalyst:')
phi = giveporosity(fitds,sigma_hc_dyn)
passd = find_d_for_first_positive(phi,fitds)
print('hc, dyn, par. passd:',passd)
#
phi = giveporosity(fitds,sigma_hc_stat)
passd = find_d_for_first_positive(phi,fitds)
print('hc, stat, par. passd:',passd)
#
phi = giveporosity(fitds,sigma_hc_dyn_perp)
passd = find_d_for_first_positive(phi,fitds)
print('hc, dyn, perp. passd:',passd)
#
phi = giveporosity(fitds,sigma_hc_stat_perp)
passd = find_d_for_first_positive(phi,fitds)
print('hc, stat, perp. passd:',passd)
print('-----------------------')
print('Cation-exchange resin membrane:')
phi = giveporosity(fitds,sigma_rm_dyn)
passd = find_d_for_first_positive(phi,fitds)
print('rm, dyn, par. passd:',passd)
#
phi = giveporosity(fitds,sigma_rm_stat)
passd = find_d_for_first_positive(phi,fitds)
print('rm, stat, par. passd:',passd)
#
phi = giveporosity(fitds,sigma_rm_dyn_perp)
passd = find_d_for_first_positive(phi,fitds)
print('rm, dyn, perp. passd:',passd)
#
phi = giveporosity(fitds,sigma_rm_stat_perp)
passd = find_d_for_first_positive(phi,fitds)
print('rm, stat,perp. passd:',passd)
print('-----------------------')

print('fitds:',fitds)

print('k_m4_stat_perp:',k_m4_stat_perp)


print('sigma_rm_stat_perp:',sigma_rm_stat_perp)
print('sigma_rm_dyn_perp:',sigma_rm_dyn_perp)
print('sigma_rm_stat:',sigma_rm_stat)
print('sigma_rm_dyn:',sigma_rm_dyn)

if showplots==True:
    plt.show()
