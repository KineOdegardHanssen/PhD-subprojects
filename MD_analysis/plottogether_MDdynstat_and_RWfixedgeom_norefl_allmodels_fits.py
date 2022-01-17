import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob

set_sigma_errors = False #True
data_sigma_errors = True # False

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

def mymodel1(d,k):
    # Assumes probability of neighbor being occupied when you're occupied: (1-k)*phi
    phi = 1 - np.pi*(1./d) # sigma fixed to 1
    teller = 1-phi*k
    nevner = phi*(1-k)
    tau = (teller/nevner)**2
    return 1./tau

def mymodel2(d,k):
    # Assumes probability of neighbor being occupied when you're occupied: vp*k
    phi = 1 - np.pi*(1./d) # sigma fixed to 1
    teller = 1+(1-phi)*(1-k)
    nevner = 1-(1-phi)*k
    tau = (teller/nevner)**2
    return 1./tau

def mymodel3(d,k):
    # Assumes probability of neighbor being occupied when you're occupied: k
    phi = 1 - np.pi*(1./d) # sigma fixed to 1
    teller = 1-phi
    nevner = 1-k
    tau = (1+teller/nevner)**2
    return 1./tau

def mymodel4(d,k,f):
    # Assumes probability of neighbor being occupied when you're occupied: kvp**f
    phi = 1 - np.pi*(1./d) # sigma fixed to 1
    vp = 1-phi
    nevner = 1-k*vp**f
    tau = (1+vp/nevner)**2
    return 1./tau

def mymodel5(d,k,f):
    # Assumes probability of neighbor being occupied when you're occupied: k+(k-1)vp**f
    phi = 1 - np.pi*(1./d) # sigma fixed to 1
    vp = 1-phi
    tau = (1+vp/(1-(k+((k-1)*vp**f))))**2
    return 1./tau

def mymodel6(d,k):
    # Assumes probability of neighbor being occupied when you're occupied: kvp**d
    phi = 1 - np.pi*(1./d) # sigma fixed to 1
    vp = 1-phi
    nevner = 1-k*vp**d
    tau = (1+vp/nevner)**2
    return 1./tau

def mymodel7(d,k):
    # Assumes probability of neighbor being occupied when you're occupied: kvp**(-1./d)
    phi = 1 - np.pi*(1./d) 
    vp = 1-phi
    nevner = 1-k*vp**(-1./d)
    tau = (1+vp/nevner)**2
    return 1./tau

def mymodel8(d,k):
    # Assumes probability of neighbor being occupied when you're occupied: k+(k-1)vp
    phi = 1 - np.pi*(1./d) # sigma fixed to 1
    vp = 1-phi
    tau = (1+vp/(1-(k+((k-1)*vp))))**2
    return 1./tau

def powerlaw1(d,k,l,n):
    return 1 - (d/l)**n*k

def powerlaw2(d,l,n):
    return 1 - (d/l)**n

def powerlaw3(d,n,k):
    return 1 - d**n*k

def powerlaw4(d,n):
    return 1 - d**n


def powerlaw3_sigmas(d,sigma,n,k):
    return 1 - (d/sigma)**n*k


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
fitds         = np.linspace(2,10,100)

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
    plotname_d_dyn     = endlocation_static+'D_vs_d_dyn_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_fittomid.png'
    plotname_d_stat    = endlocation_static+'D_vs_d_stat_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_fittomid.png'
    plotname_d_f       = endlocation_static+'D_vs_d_f_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_fittomid.png'
    # perp:
    plotname_d_dyn_perp     = endlocation_static+'D_vs_d_dyn_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_perp_fittomid.png'
    plotname_d_stat_perp    = endlocation_static+'D_vs_d_stat_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_perp_fittomid.png'
    plotname_d_f_perp       = endlocation_static+'D_vs_d_f_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_perp_fittomid.png'
    plotname_cut_d = endlocation_static+'D_vs_d_dyn_vs_stat_vs_RWgeom_cut_norefl_rsphere'+str(rsphere)+'_modelr_findsigma_fittomid.png' ## Use this?
else:
    ### d 
    plotname_d_dyn     = endlocation_static+'D_vs_d_dyn_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_fittomid.png'
    plotname_d_stat    = endlocation_static+'D_vs_d_stat_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_fittomid.png'
    plotname_d_f       = endlocation_static+'D_vs_d_f_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_fittomid.png'
    # perp:
    plotname_d_dyn_perp     = endlocation_static+'D_vs_d_dyn_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_perp_fittomid.png'
    plotname_d_stat_perp    = endlocation_static+'D_vs_d_stat_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_perp_fittomid.png'
    plotname_d_f_perp       = endlocation_static+'D_vs_d_f_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_perp_fittomid.png'
    plotname_cut_d = endlocation_static+'D_vs_d_dyn_vs_stat_cut_vs_RWgeom_big_norefl_rsphere'+str(rsphere)+'_allmodels_modelr_findsigma_fittomid.png' ## Use this?

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

indices_dyn = []
for i in range(N_dyn):
    if spacings_dyn[i]>=2 and spacings_dyn[i]<11:
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
    if spacings_stat[i]>=2 and spacings_stat[i]<11:
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
    if spacings_f[i]>=2 and spacings_f[i]<11:
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


#### Parameter estimates

cfsigmas = np.zeros(len(Dparallel_stat))+0.1
cfsigmas_f = np.zeros(len(spacings_f))+0.1
pardyn_sigmas   = Dparallel_stdv_dyn
perpdyn_sigmas  = Dzs_stdv_dyn
parstat_sigmas  = Dparallel_stdv_stat
perpstat_sigmas = Dzs_stdv_stat
# Ordered packing
#if set_sigma_errors==True:
#    popt, pcov           = curve_fit(ordered_packings, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
#if data_sigma_errors==True:
#    popt, pcov           = curve_fit(ordered_packings, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
#sigma_op_rw          = popt[0]
if set_sigma_errors==True:
    popt, pcov           = curve_fit(ordered_packings, spacings_stat, Dparallel_stat, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov           = curve_fit(ordered_packings, spacings_stat, Dparallel_stat, maxfev=thismaxfev, sigma=parstat_sigmas, absolute_sigma=True)
sigma_op_stat        = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_op_stat_rms    = rms_params[0]
if set_sigma_errors==True:
    popt, pcov           = curve_fit(ordered_packings, spacings_dyn, Dparallel_dyn, p0=7, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov           = curve_fit(ordered_packings, spacings_dyn, Dparallel_dyn, p0=7, maxfev=thismaxfev, sigma=pardyn_sigmas, absolute_sigma=True)
sigma_op_dyn         = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_op_dyn_rms     = rms_params[0]
if set_sigma_errors==True:
    popt, pcov           = curve_fit(ordered_packings, spacings_f, Dparallel_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov           = curve_fit(ordered_packings, spacings_f, Dparallel_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)
sigma_op_f             = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_op_f_rms         = rms_params[0]
# Perp
if set_sigma_errors==True:
    popt, pcov             = curve_fit(ordered_packings, spacings_stat, Dzs_stat, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov             = curve_fit(ordered_packings, spacings_stat, Dzs_stat, maxfev=thismaxfev, sigma=perpstat_sigmas, absolute_sigma=True)
sigma_op_stat_perp     = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_op_stat_perp_rms = rms_params[0]
if set_sigma_errors==True:
    popt, pcov             = curve_fit(ordered_packings, spacings_dyn, Dzs_dyn, p0=7, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov             = curve_fit(ordered_packings, spacings_dyn, Dzs_dyn, p0=7, maxfev=thismaxfev, sigma=perpdyn_sigmas, absolute_sigma=True)
sigma_op_dyn_perp      = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_op_dyn_perp_rms  = rms_params[0]
if set_sigma_errors==True:
    popt, pcov             = curve_fit(ordered_packings, spacings_f, Dzs_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov             = curve_fit(ordered_packings, spacings_f, Dzs_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)
sigma_op_f_perp        = popt[0]
rms_params             = np.sqrt(np.diag(pcov))

#sigma_op_dyn = 7
# Ds
#Ds_op_rw   = ordered_packings(fitds,sigma_op_rw)
Ds_op_stat = ordered_packings(fitds,sigma_op_stat)
Ds_op_dyn  = ordered_packings(fitds,sigma_op_dyn)
Ds_op_f    = ordered_packings(fitds,sigma_op_f)
Ds_op_stat_perp = ordered_packings(fitds,sigma_op_stat_perp)
Ds_op_dyn_perp  = ordered_packings(fitds,sigma_op_dyn_perp)
Ds_op_f_perp    = ordered_packings(fitds,sigma_op_f_perp)

# Hyperbola of revolution
#if set_sigma_errors==True:
#    popt, pcov           = curve_fit(hyp_rev, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
#if data_sigma_errors==True:
#    popt, pcov           = curve_fit(hyp_rev, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
#sigma_hr_rw          = popt[0]
if set_sigma_errors==True:
    popt, pcov           = curve_fit(hyp_rev, spacings_stat, Dparallel_stat, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov           = curve_fit(hyp_rev, spacings_stat, Dparallel_stat, maxfev=thismaxfev, sigma=parstat_sigmas, absolute_sigma=True)
sigma_hr_stat        = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_hr_stat_rms    = rms_params[0]
if set_sigma_errors==True:
    popt, pcov           = curve_fit(hyp_rev, spacings_dyn, Dparallel_dyn, p0=5, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov           = curve_fit(hyp_rev, spacings_dyn, Dparallel_dyn, p0=5, maxfev=thismaxfev, sigma=pardyn_sigmas, absolute_sigma=True)
sigma_hr_dyn           = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_hr_dyn_rms       = rms_params[0]
if set_sigma_errors==True:
    popt, pcov             = curve_fit(hyp_rev, spacings_f, Dparallel_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov             = curve_fit(hyp_rev, spacings_f, Dparallel_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)
sigma_hr_f             = popt[0]
rms_params = np.sqrt(np.diag(pcov))
# Perp
if set_sigma_errors==True:
    popt, pcov             = curve_fit(hyp_rev, spacings_stat, Dzs_stat, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov             = curve_fit(hyp_rev, spacings_stat, Dzs_stat, maxfev=thismaxfev, sigma=perpstat_sigmas, absolute_sigma=True)
sigma_hr_stat_perp     = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_hr_stat_perp_rms = rms_params[0]
if set_sigma_errors==True:
    popt, pcov             = curve_fit(hyp_rev, spacings_dyn, Dzs_dyn, p0=5, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov             = curve_fit(hyp_rev, spacings_dyn, Dzs_dyn, p0=5, maxfev=thismaxfev, sigma=perpdyn_sigmas, absolute_sigma=True)
sigma_hr_dyn_perp      = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_hr_dyn_perp_rms  = rms_params[0]
popt, pcov             = curve_fit(hyp_rev, spacings_f, Dzs_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)
sigma_hr_f_perp        = popt[0]

#sigma_hr_dyn = 5
# Ds
#Ds_hr_rw   = hyp_rev(fitds,sigma_hr_rw)
Ds_hr_stat = hyp_rev(fitds,sigma_hr_stat)
Ds_hr_dyn  = hyp_rev(fitds,sigma_hr_dyn)
Ds_hr_f    = hyp_rev(fitds,sigma_hr_f)
Ds_hr_stat_perp = hyp_rev(fitds,sigma_hr_stat_perp)
Ds_hr_dyn_perp  = hyp_rev(fitds,sigma_hr_dyn_perp)
Ds_hr_f_perp    = hyp_rev(fitds,sigma_hr_f_perp)

# Not monosized spheres
'''
popt, pcov           = curve_fit(notmonosph, spacings_rw, DRW, maxfev=thismaxfev, absolute_sigma=True)
sigma_nm_rw          = popt[0]
popt, pcov           = curve_fit(notmonosph, spacings_stat, Dparallel_stat, p0=0.5, maxfev=thismaxfev, absolute_sigma=True)
sigma_nm_stat        = popt[0]
rms_params = np.sqrt(np.diag(pcov))
popt, pcov           = curve_fit(notmonosph, spacings_dyn, Dparallel_dyn, p0=0.5, maxfev=thismaxfev, absolute_sigma=True)
sigma_nm_dyn         = popt[0]
rms_params = np.sqrt(np.diag(pcov))
popt, pcov           = curve_fit(notmonosph, spacings_f, Dparallel_f, p0=0.5, maxfev=thismaxfev, absolute_sigma=True)
sigma_nm_f           = popt[0]
# Ds
Ds_nm_rw   = notmonosph(fitds,sigma_nm_rw)
Ds_nm_stat = notmonosph(fitds,sigma_nm_stat)
Ds_nm_dyn  = notmonosph(fitds,sigma_nm_dyn)
Ds_nm_f    = notmonosph(fitds,sigma_nm_f)
# pshimd_spherepacking
popt, pcov           = curve_fit(pshimd_spherepacking, spacings_rw, DRW, maxfev=thismaxfev, absolute_sigma=True)
sigma_ps_rw          = popt[0]
popt, pcov           = curve_fit(pshimd_spherepacking, spacings_stat, Dparallel_stat, p0=0.5, maxfev=thismaxfev, absolute_sigma=True)
sigma_ps_stat        = popt[0]
rms_params = np.sqrt(np.diag(pcov))
popt, pcov           = curve_fit(pshimd_spherepacking, spacings_dyn, Dparallel_dyn, p0=0.5, maxfev=thismaxfev, absolute_sigma=True)
sigma_ps_dyn         = popt[0]
rms_params = np.sqrt(np.diag(pcov))
popt, pcov           = curve_fit(pshimd_spherepacking, spacings_f, Dparallel_f, p0=0.5, maxfev=thismaxfev, absolute_sigma=True)
sigma_ps_f           = popt[0]

# Ds
Ds_ps_rw   = pshimd_spherepacking(fitds,sigma_ps_rw)
Ds_ps_stat = pshimd_spherepacking(fitds,sigma_ps_stat)
Ds_ps_dyn  = pshimd_spherepacking(fitds,sigma_ps_dyn)
Ds_ps_f    = pshimd_spherepacking(fitds,sigma_ps_f)
# Overlapping spheres
popt, pcov           = curve_fit(overlapping_spheres, spacings_rw, DRW, maxfev=thismaxfev, absolute_sigma=True)
sigma_os_rw          = popt[0]
popt, pcov           = curve_fit(overlapping_spheres, spacings_stat, Dparallel_stat, p0=0.5, maxfev=thismaxfev, absolute_sigma=True)
sigma_os_stat        = popt[0]
rms_params = np.sqrt(np.diag(pcov))
popt, pcov           = curve_fit(overlapping_spheres, spacings_dyn, Dparallel_dyn,p0=0.5, maxfev=thismaxfev, absolute_sigma=True)
sigma_os_dyn         = popt[0]
rms_params = np.sqrt(np.diag(pcov))
popt, pcov           = curve_fit(overlapping_spheres, spacings_f, Dparallel_f, p0=0.5, maxfev=thismaxfev, absolute_sigma=True)
sigma_os_f           = popt[0]
# Ds
Ds_os_rw   = overlapping_spheres(fitds,sigma_os_rw)
Ds_os_stat = overlapping_spheres(fitds,sigma_os_stat)
Ds_os_dyn  = overlapping_spheres(fitds,sigma_os_dyn)
Ds_os_f    = overlapping_spheres(fitds,sigma_os_f)
# Overlapping cylinders
popt, pcov           = curve_fit(overlapping_cylinders, spacings_rw, DRW, maxfev=thismaxfev, absolute_sigma=True)
sigma_oc_rw          = popt[0]
popt, pcov           = curve_fit(overlapping_cylinders, spacings_stat, Dparallel_stat, sigma=sigmas_stat, p0=0.5, maxfev=thismaxfev, absolute_sigma=True)
sigma_oc_stat        = popt[0]
rms_params = np.sqrt(np.diag(pcov))
popt, pcov           = curve_fit(overlapping_cylinders, spacings_dyn, Dparallel_dyn, p0=0.5, maxfev=thismaxfev, absolute_sigma=True)
sigma_oc_dyn         = popt[0]
rms_params = np.sqrt(np.diag(pcov))
popt, pcov           = curve_fit(overlapping_cylinders, spacings_f, Dparallel_f, p0=0.5, maxfev=thismaxfev, absolute_sigma=True)
sigma_oc_f           = popt[0]
# Ds
Ds_oc_rw   = overlapping_cylinders(fitds,sigma_oc_rw)
Ds_oc_stat = overlapping_cylinders(fitds,sigma_oc_stat)
Ds_oc_dyn  = overlapping_cylinders(fitds,sigma_oc_dyn)
Ds_oc_f    = overlapping_cylinders(fitds,sigma_oc_f)
'''

# Heterogeneous catalyst
#if set_sigma_errors==True:
#    popt, pcov           = curve_fit(het_cat, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
#if data_sigma_errors==True:
#    popt, pcov           = curve_fit(het_cat, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
#sigma_hc_rw          = popt[0]
if set_sigma_errors==True:
    popt, pcov           = curve_fit(het_cat, spacings_stat, Dparallel_stat,maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov           = curve_fit(het_cat, spacings_stat, Dparallel_stat,maxfev=thismaxfev, sigma=parstat_sigmas, absolute_sigma=True)
sigma_hc_stat        = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_hc_stat_rms    = rms_params[0]
if set_sigma_errors==True:
    popt, pcov           = curve_fit(het_cat, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov           = curve_fit(het_cat, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, sigma=pardyn_sigmas, absolute_sigma=True)
sigma_hc_dyn         = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_hc_dyn_rms     = rms_params[0]
if set_sigma_errors==True:
    popt, pcov           = curve_fit(het_cat, spacings_f, Dparallel_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov           = curve_fit(het_cat, spacings_f, Dparallel_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)
sigma_hc_f           = popt[0]
# Perp
if set_sigma_errors==True:
    popt, pcov           = curve_fit(het_cat, spacings_stat, Dzs_stat,maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov           = curve_fit(het_cat, spacings_stat, Dzs_stat,maxfev=thismaxfev, sigma=perpstat_sigmas, absolute_sigma=True)
sigma_hc_stat_perp   = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_hc_stat_perp_rms = rms_params[0]
if set_sigma_errors==True:
    popt, pcov             = curve_fit(het_cat, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov             = curve_fit(het_cat, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, sigma=perpdyn_sigmas, absolute_sigma=True)
sigma_hc_dyn_perp      = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_hc_dyn_perp_rms  = rms_params[0]
if set_sigma_errors==True:
    popt, pcov             = curve_fit(het_cat, spacings_f, Dzs_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov             = curve_fit(het_cat, spacings_f, Dzs_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)
sigma_hc_f_perp        = popt[0]

# Ds
#Ds_hc_rw   = het_cat(fitds,sigma_hc_rw)
Ds_hc_stat = het_cat(fitds,sigma_hc_stat)
Ds_hc_dyn  = het_cat(fitds,sigma_hc_dyn)
Ds_hc_f    = het_cat(fitds,sigma_hc_f)
Ds_hc_stat_perp = het_cat(fitds,sigma_hc_stat_perp)
Ds_hc_dyn_perp  = het_cat(fitds,sigma_hc_dyn_perp)
Ds_hc_f_perp    = het_cat(fitds,sigma_hc_f_perp)

# Cation exchange resin membrane
#if set_sigma_errors==True:
#    popt, pcov           = curve_fit(catex_resin_membr, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
#if data_sigma_errors==True:
#    popt, pcov           = curve_fit(catex_resin_membr, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
#sigma_rm_rw          = popt[0]
if set_sigma_errors==True:
    popt, pcov           = curve_fit(catex_resin_membr, spacings_stat, Dparallel_stat, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov           = curve_fit(catex_resin_membr, spacings_stat, Dparallel_stat, maxfev=thismaxfev, sigma=parstat_sigmas, absolute_sigma=True)
sigma_rm_stat        = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_rm_stat_rms    = rms_params[0]
if set_sigma_errors==True:
    popt, pcov           = curve_fit(catex_resin_membr, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov           = curve_fit(catex_resin_membr, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, sigma=pardyn_sigmas, absolute_sigma=True)
sigma_rm_dyn         = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_rm_dyn_rms     = rms_params[0]
if set_sigma_errors==True:
    popt, pcov           = curve_fit(catex_resin_membr, spacings_f, Dparallel_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov           = curve_fit(catex_resin_membr, spacings_f, Dparallel_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)
sigma_rm_f           = popt[0]
#Perp
if set_sigma_errors==True:
    popt, pcov             = curve_fit(catex_resin_membr, spacings_stat, Dzs_stat, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov             = curve_fit(catex_resin_membr, spacings_stat, Dzs_stat, maxfev=thismaxfev, sigma=perpstat_sigmas, absolute_sigma=True)
sigma_rm_stat_perp     = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_rm_stat_perp_rms = rms_params[0]
if set_sigma_errors==True:
    popt, pcov             = curve_fit(catex_resin_membr, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov             = curve_fit(catex_resin_membr, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, sigma=perpdyn_sigmas, absolute_sigma=True)
sigma_rm_dyn_perp      = popt[0]
rms_params = np.sqrt(np.diag(pcov))
sigma_rm_dyn_perp_rms  = rms_params[0]
if set_sigma_errors==True:
    popt, pcov             = curve_fit(catex_resin_membr, spacings_f, Dzs_f, maxfev=thismaxfev, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov             = curve_fit(catex_resin_membr, spacings_f, Dzs_f, maxfev=thismaxfev, absolute_sigma=True)
sigma_rm_f_perp        = popt[0]

# Ds
#Ds_rm_rw   = catex_resin_membr(fitds,sigma_rm_rw)
Ds_rm_stat = catex_resin_membr(fitds,sigma_rm_stat)
Ds_rm_dyn  = catex_resin_membr(fitds,sigma_rm_dyn)
Ds_rm_f    = catex_resin_membr(fitds,sigma_rm_f)
Ds_rm_stat_perp = catex_resin_membr(fitds,sigma_rm_stat_perp)
Ds_rm_dyn_perp  = catex_resin_membr(fitds,sigma_rm_dyn_perp)
Ds_rm_f_perp    = catex_resin_membr(fitds,sigma_rm_f_perp)

# My models:
### mymodel1(d,k)
#if set_sigma_errors==True:
#    popt, pcov       = curve_fit(mymodel1, spacings_rw, DRW, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
#if data_sigma_errors==True:
#    popt, pcov       = curve_fit(mymodel1, spacings_rw, DRW, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
#k_m1_rw          = popt[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel1, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel1, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=(0,1), sigma=pardyn_sigmas, absolute_sigma=True)
k_m1_stat        = popt[0]
rms_params       = np.sqrt(np.diag(pcov))
k_m1_stat_stdv   = rms_params[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel1, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel1, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=(0,1), sigma=pardyn_sigmas, absolute_sigma=True)
k_m1_dyn         = popt[0]
rms_params       = np.sqrt(np.diag(pcov))
k_m1_dyn_stdv    = rms_params[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel1, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas_f)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel1, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas_f)
k_m1_f           = popt[0]
#Perp
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel1, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel1, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=(0,1), sigma=perpstat_sigmas, absolute_sigma=True)
k_m1_stat_perp      = popt[0]
rms_params          = np.sqrt(np.diag(pcov))
k_m1_stat_perp_stdv = rms_params[0]
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel1, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel1, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=(0,1), sigma=perpdyn_sigmas, absolute_sigma=True)
k_m1_dyn_perp       = popt[0]
rms_params          = np.sqrt(np.diag(pcov))
k_m1_dyn_perp_stdv  = rms_params[0]
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel1, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel1, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas_f, absolute_sigma=True)
k_m1_f_perp         = popt[0]

# Ds
#Ds_m1_rw   = mymodel1(fitds,k_m1_rw)
Ds_m1_stat = mymodel1(fitds,k_m1_stat)
Ds_m1_dyn  = mymodel1(fitds,k_m1_dyn)
Ds_m1_f    = mymodel1(fitds,k_m1_f)
Ds_m1_stat_perp = mymodel1(fitds,k_m1_stat_perp)
Ds_m1_dyn_perp  = mymodel1(fitds,k_m1_dyn_perp)
Ds_m1_f_perp    = mymodel1(fitds,k_m1_f_perp)

### mymodel2(d,k)
#if set_sigma_errors==True:
#    popt, pcov       = curve_fit(mymodel2, spacings_rw, DRW, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas, absolute_sigma=True)
#if data_sigma_errors==True:
#    popt, pcov       = curve_fit(mymodel2, spacings_rw, DRW, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas, absolute_sigma=True)
#k_m2_rw          = popt[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel2, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel2, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=(0,np.inf), sigma=parstat_sigmas, absolute_sigma=True)
k_m2_stat        = popt[0]
rms_params       = np.sqrt(np.diag(pcov))
k_m2_stat_stdv   = rms_params[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel2, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel2, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=(0,np.inf), sigma=pardyn_sigmas, absolute_sigma=True)
k_m2_dyn         = popt[0]
rms_params       = np.sqrt(np.diag(pcov))
k_m2_dyn_stdv    = rms_params[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel2, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel2, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas_f, absolute_sigma=True)
k_m2_f           = popt[0]
#Perp
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel2, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel2, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=(0,np.inf), sigma=perpstat_sigmas, absolute_sigma=True)
k_m2_stat_perp      = popt[0]
rms_params          = np.sqrt(np.diag(pcov))
k_m2_stat_perp_stdv = rms_params[0]
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel2, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel2, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=(0,np.inf), sigma=perpdyn_sigmas, absolute_sigma=True)
k_m2_dyn_perp       = popt[0]
rms_params          = np.sqrt(np.diag(pcov))
k_m2_dyn_perp_stdv  = rms_params[0]
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel2, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel2, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas_f, absolute_sigma=True)
k_m2_f_perp         = popt[0]

# Ds
#Ds_m2_rw   = mymodel2(fitds,k_m2_rw)
Ds_m2_stat = mymodel2(fitds,k_m2_stat)
Ds_m2_dyn  = mymodel2(fitds,k_m2_dyn)
Ds_m2_f    = mymodel2(fitds,k_m2_f)
Ds_m2_stat_perp = mymodel2(fitds,k_m2_stat_perp)
Ds_m2_dyn_perp  = mymodel2(fitds,k_m2_dyn_perp)
Ds_m2_f_perp    = mymodel2(fitds,k_m2_f_perp)

### mymodel3(d,k)
#if set_sigma_errors==True:
#    popt, pcov       = curve_fit(mymodel3, spacings_rw, DRW, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas, absolute_sigma=True)
#if data_sigma_errors==True:
#    popt, pcov       = curve_fit(mymodel3, spacings_rw, DRW, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas, absolute_sigma=True)
#k_m3_rw          = popt[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel3, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel3, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=(0,1), sigma=parstat_sigmas, absolute_sigma=True)
k_m3_stat        = popt[0]
rms_params       = np.sqrt(np.diag(pcov))
k_m3_stat_stdv   = rms_params[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel3, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel3, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=(0,1), sigma=pardyn_sigmas, absolute_sigma=True)
k_m3_dyn         = popt[0]
rms_params       = np.sqrt(np.diag(pcov))
k_m3_dyn_stdv    = rms_params[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel3, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel3, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas_f, absolute_sigma=True)
k_m3_f           = popt[0]
#Perp
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel3, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel3, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=(0,1), sigma=perpstat_sigmas, absolute_sigma=True)
k_m3_stat_perp      = popt[0]
rms_params          = np.sqrt(np.diag(pcov))
k_m3_stat_perp_stdv = rms_params[0]
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel3, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel3, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=(0,1), sigma=perpdyn_sigmas, absolute_sigma=True)
k_m3_dyn_perp       = popt[0]
rms_params          = np.sqrt(np.diag(pcov))
k_m3_dyn_perp_stdv  = rms_params[0]
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel3, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel3, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas_f, absolute_sigma=True)
k_m3_f_perp         = popt[0]

# Ds
#Ds_m3_rw   = mymodel3(fitds,k_m3_rw)
Ds_m3_stat = mymodel3(fitds,k_m3_stat)
Ds_m3_dyn  = mymodel3(fitds,k_m3_dyn)
Ds_m3_f    = mymodel3(fitds,k_m3_f)
Ds_m3_stat_perp = mymodel3(fitds,k_m3_stat_perp)
Ds_m3_dyn_perp  = mymodel3(fitds,k_m3_dyn_perp)
Ds_m3_f_perp    = mymodel3(fitds,k_m3_f_perp)



### mymodel4(d,k,f)
#if set_sigma_errors==True:
#    popt, pcov       = curve_fit(mymodel4, spacings_rw, DRW, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
#if data_sigma_errors==True:
#    popt, pcov       = curve_fit(mymodel4, spacings_rw, DRW, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
#k_m4_rw          = popt[0]
#f_m4_rw          = popt[1]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel4, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel4, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=(0,1), sigma=parstat_sigmas, absolute_sigma=True)
k_m4_stat        = popt[0]
f_m4_stat        = popt[1]
rms_params       = np.sqrt(np.diag(pcov))
k_m4_stat_stdv   = rms_params[0]
f_m4_stat_stdv   = rms_params[1]
pcov4_statpar    = pcov
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel4, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel4, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=(0,1), sigma=pardyn_sigmas, absolute_sigma=True)
k_m4_dyn         = popt[0]
f_m4_dyn         = popt[1]
rms_params       = np.sqrt(np.diag(pcov))
k_m4_dyn_stdv    = rms_params[0]
f_m4_dyn_stdv    = rms_params[1]
pcov4_dynpar     = pcov
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel4, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel4, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas_f, absolute_sigma=True)
k_m4_f           = popt[0]
f_m4_f           = popt[1]
#Perp
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel4, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel4, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=(0,1), sigma=perpstat_sigmas, absolute_sigma=True)
k_m4_stat_perp      = popt[0]
f_m4_stat_perp      = popt[1]
rms_params          = np.sqrt(np.diag(pcov))
k_m4_stat_perp_stdv = rms_params[0]
f_m4_stat_perp_stdv = rms_params[1]
pcov4_statperp      = pcov
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel4, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel4, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=(0,1), sigma=perpdyn_sigmas, absolute_sigma=True)
k_m4_dyn_perp       = popt[0]
f_m4_dyn_perp       = popt[1]
rms_params          = np.sqrt(np.diag(pcov))
k_m4_dyn_perp_stdv  = rms_params[0]
f_m4_dyn_perp_stdv  = rms_params[1]
pcov4_dynperp       = pcov
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel4, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel4, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas_f, absolute_sigma=True)
k_m4_f_perp         = popt[0]
f_m4_f_perp         = popt[1]

# Ds
#Ds_m4_rw   = mymodel4(fitds,k_m4_rw,f_m4_rw)
Ds_m4_stat = mymodel4(fitds,k_m4_stat,f_m4_stat)
Ds_m4_dyn  = mymodel4(fitds,k_m4_dyn,f_m4_dyn)
Ds_m4_f    = mymodel4(fitds,k_m4_f,f_m4_f)
Ds_m4_stat_perp = mymodel4(fitds,k_m4_stat_perp,f_m4_stat_perp)
Ds_m4_dyn_perp  = mymodel4(fitds,k_m4_dyn_perp,f_m4_dyn_perp)
Ds_m4_f_perp    = mymodel4(fitds,k_m4_f_perp,f_m4_f_perp)


#'''
### mymodel5(d,k,f) # Doesn't work.
#if set_sigma_errors==True:
#    popt, pcov      = curve_fit(mymodel5, spacings_rw, DRW, maxfev=1000000, sigma=cfsigmas, absolute_sigma=True)
#if data_sigma_errors==True:
#    popt, pcov      = curve_fit(mymodel5, spacings_rw, DRW, maxfev=1000000, sigma=cfsigmas, absolute_sigma=True)
#k_m5_rw         = popt[0]
#f_m5_rw         = popt[1]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel5, spacings_stat, Dparallel_stat, maxfev=1000000, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel5, spacings_stat, Dparallel_stat, maxfev=1000000, sigma=parstat_sigmas, absolute_sigma=True)
k_m5_stat        = popt[0]
f_m5_stat        = popt[1]
rms_params       = np.sqrt(np.diag(pcov))
k_m5_stat_stdv   = rms_params[0]
f_m5_stat_stdv   = rms_params[1]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel5, spacings_dyn, Dparallel_dyn, maxfev=1000000, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel5, spacings_dyn, Dparallel_dyn, maxfev=1000000, sigma=pardyn_sigmas, absolute_sigma=True)
k_m5_dyn         = popt[0]
f_m5_dyn         = popt[1]
rms_params       = np.sqrt(np.diag(pcov))
k_m5_dyn_stdv    = rms_params[0]
f_m5_dyn_stdv    = rms_params[1]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel5, spacings_f, Dparallel_f, maxfev=1000000, sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel5, spacings_f, Dparallel_f, maxfev=1000000, sigma=cfsigmas_f, absolute_sigma=True)
k_m5_f           = popt[0]
f_m5_f           = popt[1]
#Perp
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel5, spacings_stat, Dzs_stat, maxfev=1000000, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel5, spacings_stat, Dzs_stat, maxfev=1000000, sigma=perpstat_sigmas, absolute_sigma=True)
k_m5_stat_perp      = popt[0]
f_m5_stat_perp      = popt[1]
rms_params          = np.sqrt(np.diag(pcov))
k_m5_stat_perp_stdv = rms_params[0]
f_m5_stat_perp_stdv = rms_params[1]
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel5, spacings_dyn, Dzs_dyn, maxfev=1000000, sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel5, spacings_dyn, Dzs_dyn, maxfev=1000000, sigma=perpdyn_sigmas, absolute_sigma=True)
k_m5_dyn_perp       = popt[0]
f_m5_dyn_perp       = popt[1]
rms_params          = np.sqrt(np.diag(pcov))
k_m5_dyn_perp_stdv  = rms_params[0]
f_m5_dyn_perp_stdv  = rms_params[1]
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel5, spacings_f, Dzs_f, maxfev=1000000, sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel5, spacings_f, Dzs_f, maxfev=1000000, sigma=cfsigmas_f, absolute_sigma=True)
k_m5_f_perp         = popt[0]
f_m5_f_perp         = popt[1]

# Ds
#Ds_m5_rw   = mymodel5(fitds,k_m5_rw,f_m5_rw)
Ds_m5_stat = mymodel5(fitds,k_m5_stat,f_m5_stat)
Ds_m5_dyn  = mymodel5(fitds,k_m5_dyn,f_m5_dyn)
Ds_m5_f    = mymodel5(fitds,k_m5_f,f_m5_f)
Ds_m5_stat_perp = mymodel5(fitds,k_m5_stat_perp,f_m5_stat_perp)
Ds_m5_dyn_perp  = mymodel5(fitds,k_m5_dyn_perp,f_m5_dyn_perp)
Ds_m5_f_perp    = mymodel5(fitds,k_m5_f_perp,f_m5_f_perp)
#'''

### mymodel6(d,k)
#if set_sigma_errors==True:
#    popt, pcov       = curve_fit(mymodel6, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,np.inf))
#if data_sigma_errors==True:
#    popt, pcov       = curve_fit(mymodel6, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,np.inf))
#k_m6_rw          = popt[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel6, spacings_stat, Dparallel_stat, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,np.inf))
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel6, spacings_stat, Dparallel_stat, maxfev=thismaxfev, sigma=parstat_sigmas, absolute_sigma=True)#, bounds=(0,np.inf))
k_m6_stat        = popt[0]
rms_params       = np.sqrt(np.diag(pcov))
k_m6_stat_stdv   = rms_params[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel6, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,np.inf))
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel6, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, sigma=pardyn_sigmas, absolute_sigma=True)#, bounds=(0,np.inf))
k_m6_dyn         = popt[0]
rms_params       = np.sqrt(np.diag(pcov))
k_m6_dyn_stdv    = rms_params[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel6, spacings_f, Dparallel_f, sigma=cfsigmas_f, maxfev=thismaxfev)#, bounds=(0,np.inf))
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel6, spacings_f, Dparallel_f, sigma=cfsigmas_f, maxfev=thismaxfev)#, bounds=(0,np.inf))
k_m6_f           = popt[0]
#Perp
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel6, spacings_stat, Dzs_stat, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,np.inf))
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel6, spacings_stat, Dzs_stat, maxfev=thismaxfev, sigma=perpstat_sigmas, absolute_sigma=True)#, bounds=(0,np.inf))
k_m6_stat_perp      = popt[0]
rms_params          = np.sqrt(np.diag(pcov))
k_m6_stat_perp_stdv = rms_params[0]
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel6, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,np.inf))
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel6, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, sigma=perpstat_sigmas, absolute_sigma=True)#, bounds=(0,np.inf))
k_m6_dyn_perp       = popt[0]
rms_params          = np.sqrt(np.diag(pcov))
k_m6_dyn_perp_stdv  = rms_params[0]
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel6, spacings_f, Dzs_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,np.inf))
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel6, spacings_f, Dzs_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,np.inf))
k_m6_f_perp         = popt[0]

# Ds
#Ds_m6_rw   = mymodel6(fitds,k_m6_rw)
Ds_m6_stat = mymodel6(fitds,k_m6_stat)
Ds_m6_dyn  = mymodel6(fitds,k_m6_dyn)
Ds_m6_f    = mymodel6(fitds,k_m6_f)
Ds_m6_stat_perp = mymodel6(fitds,k_m6_stat_perp)
Ds_m6_dyn_perp  = mymodel6(fitds,k_m6_dyn_perp)
Ds_m6_f_perp    = mymodel6(fitds,k_m6_f_perp)

### mymodel7(d,k)
#if set_sigma_errors==True:
#    popt, pcov         = curve_fit(mymodel7, spacings_rw, DRW, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas, absolute_sigma=True)
#if data_sigma_errors==True:
#    popt, pcov         = curve_fit(mymodel7, spacings_rw, DRW, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas, absolute_sigma=True)
#k_m7_rw            = popt[0]
if set_sigma_errors==True:
    popt, pcov         = curve_fit(mymodel7, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov         = curve_fit(mymodel7, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=(0,np.inf), sigma=parstat_sigmas, absolute_sigma=True)
rms_params = np.sqrt(np.diag(pcov))
k_m7_stat          = popt[0]
k_m7_stat_stdv     = rms_params[0]
if set_sigma_errors==True:
    popt, pcov         = curve_fit(mymodel7, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov         = curve_fit(mymodel7, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=(0,np.inf), sigma=pardyn_sigmas, absolute_sigma=True)
rms_params = np.sqrt(np.diag(pcov))
k_m7_dyn           = popt[0]
k_m7_dyn_stdv      = rms_params[0]
if set_sigma_errors==True:
    popt, pcov         = curve_fit(mymodel7, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov         = curve_fit(mymodel7, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas_f, absolute_sigma=True)
k_m7_f             = popt[0]
#Perp
if set_sigma_errors==True:
    popt, pcov              = curve_fit(mymodel7, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov              = curve_fit(mymodel7, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=(0,np.inf), sigma=perpstat_sigmas, absolute_sigma=True)
rms_params = np.sqrt(np.diag(pcov))
k_m7_stat_perp          = popt[0]
k_m7_stat_perp_stdv     = rms_params[0]
if set_sigma_errors==True:
    popt, pcov              = curve_fit(mymodel7, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov              = curve_fit(mymodel7, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=(0,np.inf), sigma=perpdyn_sigmas, absolute_sigma=True)
rms_params = np.sqrt(np.diag(pcov))
k_m7_dyn_perp           = popt[0]
k_m7_dyn_perp_stdv      = rms_params[0]
if set_sigma_errors==True:
    popt, pcov              = curve_fit(mymodel7, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov              = curve_fit(mymodel7, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=(0,np.inf), sigma=cfsigmas_f, absolute_sigma=True)
k_m7_f_perp             = popt[0]

# Ds
#Ds_m7_rw   = mymodel7(fitds,k_m7_rw)
Ds_m7_stat = mymodel7(fitds,k_m7_stat)
Ds_m7_dyn  = mymodel7(fitds,k_m7_dyn)
Ds_m7_f    = mymodel7(fitds,k_m7_f)
Ds_m7_stat_perp = mymodel7(fitds,k_m7_stat_perp)
Ds_m7_dyn_perp  = mymodel7(fitds,k_m7_dyn_perp)
Ds_m7_f_perp    = mymodel7(fitds,k_m7_f_perp)

### 
### mymodel8(d,k)
#if set_sigma_errors==True:
#    popt, pcov       = curve_fit(mymodel8, spacings_rw, DRW, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
#if data_sigma_errors==True:
#    popt, pcov       = curve_fit(mymodel8, spacings_rw, DRW, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
#k_m8_rw          = popt[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel8, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel8, spacings_stat, Dparallel_stat, maxfev=thismaxfev, bounds=(0,1), sigma=parstat_sigmas, absolute_sigma=True)
k_m8_stat        = popt[0]
rms_params       = np.sqrt(np.diag(pcov))
k_m8_stat_stdv   = rms_params[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel8, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel8, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, bounds=(0,1), sigma=pardyn_sigmas, absolute_sigma=True)
k_m8_dyn         = popt[0]
rms_params       = np.sqrt(np.diag(pcov))
k_m8_dyn_stdv    = rms_params[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel8, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov       = curve_fit(mymodel8, spacings_f, Dparallel_f, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas_f, absolute_sigma=True)
k_m8_f           = popt[0]
#Perp
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel8, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel8, spacings_stat, Dzs_stat, maxfev=thismaxfev, bounds=(0,1), sigma=perpstat_sigmas, absolute_sigma=True)
k_m8_stat_perp      = popt[0]
rms_params          = np.sqrt(np.diag(pcov))
k_m8_stat_perp_stdv = rms_params[0]
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel8, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel8, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, bounds=(0,1), sigma=perpdyn_sigmas, absolute_sigma=True)
k_m8_dyn_perp       = popt[0]
rms_params          = np.sqrt(np.diag(pcov))
k_m8_dyn_perp_stdv  = rms_params[0]
if set_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel8, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas_f, absolute_sigma=True)
if data_sigma_errors==True:
    popt, pcov          = curve_fit(mymodel8, spacings_f, Dzs_f, maxfev=thismaxfev, bounds=(0,1), sigma=cfsigmas_f, absolute_sigma=True)
k_m8_f_perp         = popt[0]

# Ds
#Ds_m8_rw   = mymodel8(fitds,k_m8_rw)
Ds_m8_stat = mymodel8(fitds,k_m8_stat)
Ds_m8_dyn  = mymodel8(fitds,k_m8_dyn)
Ds_m8_f    = mymodel8(fitds,k_m8_f)
Ds_m8_stat_perp = mymodel8(fitds,k_m8_stat_perp)
Ds_m8_dyn_perp  = mymodel8(fitds,k_m8_dyn_perp)
Ds_m8_f_perp    = mymodel8(fitds,k_m8_f_perp)



### powerlaw1(d,k,l,n)
#if set_sigma_errors==True:
#    popt, pcov       = curve_fit(powerlaw1, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
#if data_sigma_errors==True:
#    popt, pcov       = curve_fit(powerlaw1, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
#k_p1_rw          = popt[0]
#l_p1_rw          = popt[1]
#n_p1_rw          = popt[2]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw1, spacings_stat, Dparallel_stat, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw1, spacings_stat, Dparallel_stat, maxfev=thismaxfev, sigma=parstat_sigmas, absolute_sigma=True)#, bounds=(0,1))
k_p1_stat        = popt[0]
l_p1_stat        = popt[1]
n_p1_stat        = popt[2]
rms_params       = np.sqrt(np.diag(pcov))
k_p1_stat_stdv   = rms_params[0]
l_p1_stat_stdv   = rms_params[1]
n_p1_stat_stdv   = rms_params[2]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw1, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw1, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, sigma=pardyn_sigmas, absolute_sigma=True)#, bounds=(0,1))
k_p1_dyn         = popt[0]
l_p1_dyn         = popt[1]
n_p1_dyn         = popt[2]
rms_params       = np.sqrt(np.diag(pcov))
k_p1_dyn_stdv    = rms_params[0]
l_p1_dyn_stdv    = rms_params[1]
n_p1_dyn_stdv    = rms_params[2]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw1, spacings_f, Dparallel_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw1, spacings_f, Dparallel_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,1))
k_p1_f           = popt[0]
l_p1_f           = popt[1]
n_p1_f           = popt[2]
#Perp
if set_sigma_errors==True:
    popt, pcov          = curve_fit(powerlaw1, spacings_stat, Dzs_stat, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov          = curve_fit(powerlaw1, spacings_stat, Dzs_stat, maxfev=thismaxfev, sigma=perpstat_sigmas, absolute_sigma=True)#, bounds=(0,1))
k_p1_stat_perp        = popt[0]
l_p1_stat_perp        = popt[1]
n_p1_stat_perp        = popt[2]
rms_params            = np.sqrt(np.diag(pcov))
k_p1_stat_perp_stdv   = rms_params[0]
l_p1_stat_perp_stdv   = rms_params[1]
n_p1_stat_perp_stdv   = rms_params[2]
if set_sigma_errors==True:
    popt, pcov            = curve_fit(powerlaw1, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov            = curve_fit(powerlaw1, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, sigma=perpdyn_sigmas, absolute_sigma=True)#, bounds=(0,1))
k_p1_dyn_perp        = popt[0]
l_p1_dyn_perp        = popt[1]
n_p1_dyn_perp        = popt[2]
rms_params           = np.sqrt(np.diag(pcov))
k_p1_dyn_perp_stdv   = rms_params[0]
l_p1_dyn_perp_stdv   = rms_params[1]
n_p1_dyn_perp_stdv   = rms_params[2]
if set_sigma_errors==True:
    popt, pcov           = curve_fit(powerlaw1, spacings_f, Dzs_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov           = curve_fit(powerlaw1, spacings_f, Dzs_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,1))
k_p1_f_perp         = popt[0]
l_p1_f_perp         = popt[1]
n_p1_f_perp         = popt[2]

# Ds
#Ds_p1_rw   = powerlaw1(fitds,k_p1_rw,l_p1_rw,n_p1_rw)
Ds_p1_stat = powerlaw1(fitds,k_p1_stat,l_p1_stat,n_p1_stat)
Ds_p1_dyn  = powerlaw1(fitds,k_p1_dyn,l_p1_dyn,n_p1_dyn)
Ds_p1_f    = powerlaw1(fitds,k_p1_f,l_p1_f,n_p1_f)
Ds_p1_stat_perp = powerlaw1(fitds,k_p1_stat_perp,l_p1_stat_perp,n_p1_stat_perp)
Ds_p1_dyn_perp  = powerlaw1(fitds,k_p1_dyn_perp,l_p1_dyn_perp,n_p1_dyn_perp)
Ds_p1_f_perp    = powerlaw1(fitds,k_p1_f_perp,l_p1_f_perp,n_p1_f_perp)

##### powerlaw2(d,l,n)
#if set_sigma_errors==True:
#    popt, pcov       = curve_fit(powerlaw2, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
#if data_sigma_errors==True:
#    popt, pcov       = curve_fit(powerlaw2, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
#l_p2_rw          = popt[0]
#n_p2_rw          = popt[1]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw2, spacings_stat, Dparallel_stat, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw2, spacings_stat, Dparallel_stat, maxfev=thismaxfev, sigma=parstat_sigmas, absolute_sigma=True)#, bounds=(0,1))
l_p2_stat        = popt[0]
n_p2_stat        = popt[1]
rms_params       = np.sqrt(np.diag(pcov))
l_p2_stat_stdv   = rms_params[0]
n_p2_stat_stdv   = rms_params[1]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw2, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw2, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, sigma=pardyn_sigmas, absolute_sigma=True)#, bounds=(0,1))
l_p2_dyn         = popt[0]
n_p2_dyn         = popt[1]
rms_params       = np.sqrt(np.diag(pcov))
l_p2_dyn_stdv    = rms_params[0]
n_p2_dyn_stdv    = rms_params[1]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw2, spacings_f, Dparallel_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw2, spacings_f, Dparallel_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,1))
l_p2_f           = popt[0]
n_p2_f           = popt[1]
#Perp
if set_sigma_errors==True:
    popt, pcov          = curve_fit(powerlaw2, spacings_stat, Dzs_stat, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov          = curve_fit(powerlaw2, spacings_stat, Dzs_stat, maxfev=thismaxfev, sigma=perpstat_sigmas, absolute_sigma=True)#, bounds=(0,1))
l_p2_stat_perp        = popt[0]
n_p2_stat_perp        = popt[1]
rms_params            = np.sqrt(np.diag(pcov))
l_p2_stat_perp_stdv   = rms_params[0]
n_p2_stat_perp_stdv   = rms_params[1]
if set_sigma_errors==True:
    popt, pcov            = curve_fit(powerlaw2, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov            = curve_fit(powerlaw2, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, sigma=perpdyn_sigmas, absolute_sigma=True)#, bounds=(0,1))
l_p2_dyn_perp        = popt[0]
n_p2_dyn_perp        = popt[1]
rms_params           = np.sqrt(np.diag(pcov))
l_p2_dyn_perp_stdv   = rms_params[0]
n_p2_dyn_perp_stdv   = rms_params[1]
if set_sigma_errors==True:
    popt, pcov           = curve_fit(powerlaw2, spacings_f, Dzs_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov           = curve_fit(powerlaw2, spacings_f, Dzs_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,1))
l_p2_f_perp         = popt[0]
n_p2_f_perp         = popt[1]

# Ds
#Ds_p2_rw   = powerlaw2(fitds,l_p2_rw,n_p2_rw)
Ds_p2_stat = powerlaw2(fitds,l_p2_stat,n_p2_stat)
Ds_p2_dyn  = powerlaw2(fitds,l_p2_dyn,n_p2_dyn)
Ds_p2_f    = powerlaw2(fitds,l_p2_f,n_p2_f)
Ds_p2_stat_perp = powerlaw2(fitds,l_p2_stat_perp,n_p2_stat_perp)
Ds_p2_dyn_perp  = powerlaw2(fitds,l_p2_dyn_perp,n_p2_dyn_perp)
Ds_p2_f_perp    = powerlaw2(fitds,l_p2_f_perp,n_p2_f_perp)

##### powerlaw3(d,n,k)
#if set_sigma_errors==True:
#    popt, pcov       = curve_fit(powerlaw3, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
#if data_sigma_errors==True:
#    popt, pcov       = curve_fit(powerlaw3, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
#n_p3_rw          = popt[0]
#k_p3_rw          = popt[1]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw3, spacings_stat, Dparallel_stat, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw3, spacings_stat, Dparallel_stat, maxfev=thismaxfev, sigma=parstat_sigmas, absolute_sigma=True)#, bounds=(0,1))
pcovp3_statpar   = pcov
n_p3_stat        = popt[0]
k_p3_stat        = popt[1]
rms_params       = np.sqrt(np.diag(pcov))
n_p3_stat_stdv   = rms_params[0]
k_p3_stat_stdv   = rms_params[1]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw3, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw3, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, sigma=pardyn_sigmas, absolute_sigma=True)#, bounds=(0,1))
pcovp3_dynpar    = pcov
n_p3_dyn         = popt[0]
k_p3_dyn         = popt[1]
rms_params       = np.sqrt(np.diag(pcov))
n_p3_dyn_stdv    = rms_params[0]
k_p3_dyn_stdv    = rms_params[1]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw3, spacings_f, Dparallel_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw3, spacings_f, Dparallel_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,1))
n_p3_f           = popt[0]
k_p3_f           = popt[1]
#Perp
if set_sigma_errors==True:
    popt, pcov            = curve_fit(powerlaw3, spacings_stat, Dzs_stat, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov            = curve_fit(powerlaw3, spacings_stat, Dzs_stat, maxfev=thismaxfev, sigma=perpstat_sigmas, absolute_sigma=True)#, bounds=(0,1))
pcovp3_statperp       = pcov
n_p3_stat_perp        = popt[0]
k_p3_stat_perp        = popt[1]
rms_params            = np.sqrt(np.diag(pcov))
n_p3_stat_perp_stdv   = rms_params[0]
k_p3_stat_perp_stdv   = rms_params[1]
if set_sigma_errors==True:
    popt, pcov            = curve_fit(powerlaw3, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov            = curve_fit(powerlaw3, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, sigma=perpdyn_sigmas, absolute_sigma=True)#, bounds=(0,1))
pcovp3_dynperp        = pcov
n_p3_dyn_perp         = popt[0]
k_p3_dyn_perp         = popt[1]
rms_params            = np.sqrt(np.diag(pcov))
n_p3_dyn_perp_stdv    = rms_params[0]
k_p3_dyn_perp_stdv    = rms_params[1]
if set_sigma_errors==True:
    popt, pcov            = curve_fit(powerlaw3, spacings_f, Dzs_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov            = curve_fit(powerlaw3, spacings_f, Dzs_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,1))
n_p3_f_perp           = popt[0]
k_p3_f_perp           = popt[1]

# Ds
#Ds_p3_rw   = powerlaw3(fitds,n_p3_rw,k_p3_rw)
Ds_p3_stat = powerlaw3(fitds,n_p3_stat,k_p3_stat)
Ds_p3_dyn  = powerlaw3(fitds,n_p3_dyn,k_p3_dyn)
Ds_p3_f    = powerlaw3(fitds,n_p3_f,k_p3_f)
Ds_p3_stat_perp = powerlaw3(fitds,n_p3_stat_perp,k_p3_stat_perp)
Ds_p3_dyn_perp  = powerlaw3(fitds,n_p3_dyn_perp,k_p3_dyn_perp)
Ds_p3_f_perp    = powerlaw3(fitds,n_p3_f_perp,k_p3_f_perp)


##### powerlaw4(d,n)
#if set_sigma_errors==True:
#    popt, pcov       = curve_fit(powerlaw4, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
#if data_sigma_errors==True:
#    popt, pcov       = curve_fit(powerlaw4, spacings_rw, DRW, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
#n_p4_rw          = popt[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw4, spacings_stat, Dparallel_stat, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw4, spacings_stat, Dparallel_stat, maxfev=thismaxfev, sigma=parstat_sigmas, absolute_sigma=True)#, bounds=(0,1))
n_p4_stat        = popt[0]
rms_params       = np.sqrt(np.diag(pcov))
n_p4_stat_stdv   = rms_params[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw4, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw4, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev, sigma=pardyn_sigmas, absolute_sigma=True)#, bounds=(0,1))
n_p4_dyn         = popt[0]
rms_params       = np.sqrt(np.diag(pcov))
n_p4_dyn_stdv    = rms_params[0]
if set_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw4, spacings_f, Dparallel_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov       = curve_fit(powerlaw4, spacings_f, Dparallel_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,1))
n_p4_f           = popt[0]
#Perp
if set_sigma_errors==True:
    popt, pcov          = curve_fit(powerlaw4, spacings_stat, Dzs_stat, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov          = curve_fit(powerlaw4, spacings_stat, Dzs_stat, maxfev=thismaxfev, sigma=perpstat_sigmas, absolute_sigma=True)#, bounds=(0,1))
n_p4_stat_perp       = popt[0]
rms_params           = np.sqrt(np.diag(pcov))
n_p4_stat_perp_stdv  = rms_params[0]
if set_sigma_errors==True:
    popt, pcov           = curve_fit(powerlaw4, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, sigma=cfsigmas, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov           = curve_fit(powerlaw4, spacings_dyn, Dzs_dyn, maxfev=thismaxfev, sigma=perpdyn_sigmas, absolute_sigma=True)#, bounds=(0,1))
n_p4_dyn_perp        = popt[0]
rms_params           = np.sqrt(np.diag(pcov))
n_p4_dyn_perp_stdv   = rms_params[0]
if set_sigma_errors==True:
    popt, pcov           = curve_fit(powerlaw4, spacings_f, Dzs_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,1))
if data_sigma_errors==True:
    popt, pcov           = curve_fit(powerlaw4, spacings_f, Dzs_f, maxfev=thismaxfev, sigma=cfsigmas_f, absolute_sigma=True)#, bounds=(0,1))
n_p4_f_perp         = popt[0]

# Ds
#Ds_p4_rw   = powerlaw4(fitds,n_p4_rw)
Ds_p4_stat = powerlaw4(fitds,n_p4_stat)
Ds_p4_dyn  = powerlaw4(fitds,n_p4_dyn)
Ds_p4_f    = powerlaw4(fitds,n_p4_f)
Ds_p4_stat_perp = powerlaw4(fitds,n_p4_stat_perp)
Ds_p4_dyn_perp  = powerlaw4(fitds,n_p4_dyn_perp)
Ds_p4_f_perp    = powerlaw4(fitds,n_p4_f_perp)



#### Plotting:
##params = {'mathtext.default': 'regular' }  # Ditch this.   
###plt.rcParams.update(params)
#'''
plt.figure(figsize=(6,5))#(14,5))
ax = plt.subplot(111)
ax.plot(spacings_dyn, Dparallel_dyn, color='g', label=r'$D_\parallel$, dyn.')
ax.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
#ax.plot(dg,DRW, color='b', label='Random walk, fixed geom., no refl.')
#ax.plot(fitds,Ds_op_dyn, '--', label=r'Ordered packings') # I have most values for the dynamic brush
#ax.plot(fitds,Ds_hr_dyn, '-.', label=r'Hyperbola of revolution') # 2
#ax.plot(fitds,Ds_nm_dyn, '--', label=r'Not monosized spheres')  ####
#ax.plot(fitds,Ds_ps_dyn, '--', label=r'Monodisperse sphere packing') ####
#ax.plot(fitds,Ds_os_dyn, '--', label=r'Overlapping spheres') # 5 ####
#ax.plot(fitds,Ds_oc_dyn, '--', label=r'Overlapping cylinders')  ####
#ax.plot(fitds,Ds_hc_dyn, ':', label=r'Heterogeneous catalyst') 
#ax.plot(fitds,Ds_m1_dyn, '-.', label=r'Custom model 1')
#ax.plot(fitds,Ds_m2_dyn, '-.', label=r'Custom model 2')
#ax.plot(fitds,Ds_m3_dyn, '-.', label=r'Custom model 3')
ax.plot(fitds,Ds_m4_dyn, '-.', label=r'Custom model')
#ax.plot(fitds,Ds_m5_dyn, '-.', label=r'Custom model 5')
#ax.plot(fitds,Ds_m6_dyn, '-.', label=r'Custom model 6')
#ax.plot(fitds,Ds_m7_dyn, '-.', label=r'Custom model 7')
#ax.plot(fitds,Ds_m8_dyn, '-.', label=r'Custom model 8')
#ax.plot(fitds,Ds_p1_dyn, '-.', label=r'Power law 1')
#ax.plot(fitds,Ds_p2_dyn, '-.', label=r'Power law 2')
ax.plot(fitds,Ds_p3_dyn, '--', label=r'Power law')
#ax.plot(fitds,Ds_p4_dyn, '-.', label=r'Power law 4')
line1, = ax.plot(fitds,Ds_rm_dyn, '--,', label=r'Cation-exchange resin membrane') 
line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
# Add more plotting! # Need to beware of negative porosity for some!
plt.xlabel(r'$d$ (nm)', fontsize=14)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=14)
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.legend(loc='upper left', fontsize=11)
plt.tight_layout()
plt.savefig(plotname_d_dyn)


'''
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
#ax.plot(fitds,Ds_hc_stat, ':', label=r'Heterogeneous catalyst') 
#ax.plot(fitds,Ds_m1_stat, '-.', label=r'Custom model 1')
#ax.plot(fitds,Ds_m2_stat, '-.', label=r'Custom model 2')
#ax.plot(fitds,Ds_m3_stat, '-.', label=r'Custom model 3')
ax.plot(fitds,Ds_m4_stat, '-.', label=r'Custom model')
#ax.plot(fitds,Ds_m5_stat, '-.', label=r'Custom model 5')
#ax.plot(fitds,Ds_m6_stat, '-.', label=r'Custom model 6')
#ax.plot(fitds,Ds_m7_stat, '-.', label=r'Custom model 7')
#ax.plot(fitds,Ds_m8_stat, '-.', label=r'Custom model 8')
#ax.plot(fitds,Ds_p1_stat, '-.', label=r'Power law 1')
#ax.plot(fitds,Ds_p2_stat, '-.', label=r'Power law 2')
ax.plot(fitds,Ds_p3_stat, '--', label=r'Power law')
#ax.plot(fitds,Ds_p4_stat, '-.', label=r'Power law 4')
line1, = ax.plot(fitds,Ds_rm_stat, '--,', label=r'Cation-exchange resin membrane') 
line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
# Add more plotting! # Need to beware of negative porosity for some!
plt.xlabel(r'$d$ (nm)', fontsize=14)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=14)
plt.legend(loc='lower right')
plt.tight_layout()
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_d_stat)
'''

#... Do the same, men for perp.
plt.figure(figsize=(6,5))#(14,5))
ax = plt.subplot(111)
ax.plot(spacings_dyn, Dzs_dyn, color='g', label=r'$D_\perp$, dyn.')
ax.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='g', alpha=0.2)
#ax.plot(fitds,Ds_op_dyn_perp, '--', label=r'Ordered packings') # I have most values for the dynamic brush
#ax.plot(fitds,Ds_hr_dyn_perp, '-.', label=r'Hyperbola of revolution') # 2
#ax.plot(fitds,Ds_nm_dyn_perp, '--', label=r'Not monosized spheres')  ####
#ax.plot(fitds,Ds_ps_dyn_perp, '--', label=r'Monodisperse sphere packing') ####
#ax.plot(fitds,Ds_os_dyn_perp, '--', label=r'Overlapping spheres') # 5 ####
#ax.plot(fitds,Ds_oc_dyn_perp, '--', label=r'Overlapping cylinders')  ####
#ax.plot(fitds,Ds_hc_dyn_perp, ':', label=r'Heterogeneous catalyst') 
#ax.plot(fitds,Ds_m1_dyn_perp, '-.', label=r'Custom model 1')
#ax.plot(fitds,Ds_m2_dyn_perp, '-.', label=r'Custom model 2')
#ax.plot(fitds,Ds_m3_dyn_perp, '-.', label=r'Custom model 3')
ax.plot(fitds,Ds_m4_dyn_perp, '-.', label=r'Custom model') #4')
#ax.plot(fitds,Ds_m5_dyn_perp, '-.', label=r'Custom model 5')
#ax.plot(fitds,Ds_m6_dyn_perp, '-.', label=r'Custom model 6')
#ax.plot(fitds,Ds_m7_dyn_perp, '-.', label=r'Custom model 7')
#ax.plot(fitds,Ds_m8_dyn_perp, '-.', label=r'Custom model 8')
#ax.plot(fitds,Ds_p1_dyn_perp, '-.', label=r'Power law 1')
#ax.plot(fitds,Ds_p2_dyn_perp, '-.', label=r'Power law 2')
ax.plot(fitds,Ds_p3_dyn_perp, '--', label=r'Power law')# 3')
#ax.plot(fitds,Ds_p4_dyn_perp, '-.', label=r'Power law 4')
line1, = ax.plot(fitds,Ds_rm_dyn_perp, '--,', label=r'Cation-exchange resin membrane') 
line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
# Add more plotting! # Need to beware of negative porosity for some!
plt.xlabel(r'$d$ (nm)', fontsize=14)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=14)
plt.legend(loc='lower right', fontsize=11)
plt.tight_layout()
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_d_dyn_perp)

'''
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
#ax.plot(fitds,Ds_hc_stat_perp, ':', label=r'Heterogeneous catalyst') 
#ax.plot(fitds,Ds_m1_stat_perp, '-.', label=r'Custom model 1')
#ax.plot(fitds,Ds_m2_stat_perp, '-.', label=r'Custom model 2')
#ax.plot(fitds,Ds_m3_stat_perp, '-.', label=r'Custom model 3')
ax.plot(fitds,Ds_m4_stat_perp, '-.', label=r'Custom model') #4')
#ax.plot(fitds,Ds_m5_stat_perp, '-.', label=r'Custom model 5')
#ax.plot(fitds,Ds_m6_stat_perp, '-.', label=r'Custom model 6')
#ax.plot(fitds,Ds_m7_stat_perp, '-.', label=r'Custom model 7')
#ax.plot(fitds,Ds_m8_stat_perp, '-.', label=r'Custom model 8')
#ax.plot(fitds,Ds_p1_stat_perp, '-.', label=r'Power law 1')
#ax.plot(fitds,Ds_p2_stat_perp, '-.', label=r'Power law 2')
ax.plot(fitds,Ds_p3_stat_perp, '--', label=r'Power law') # 3')
#ax.plot(fitds,Ds_p4_stat_perp, '-.', label=r'Power law 4')
line1, = ax.plot(fitds,Ds_rm_stat_perp, '--,', label=r'Cation-exchange resin membrane') 
line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
# Add more plotting! # Need to beware of negative porosity for some!
plt.xlabel(r'$d$ (nm)', fontsize=14)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=14)
plt.legend(loc='lower right')
plt.tight_layout()
#plt.title(r'$D/D_{\mathregular{bulk}}$ vs $d$')
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
#plt.tight_layout()
plt.savefig(plotname_d_stat_perp)

#### Plotting +/- stdv: #####

## Model 2
# Parallel
Ds_m2_stat_p = mymodel2(fitds,k_m2_stat+k_m2_stat_stdv)
Ds_m2_stat_m = mymodel2(fitds,k_m2_stat-k_m2_stat_stdv)
Ds_m2_dyn_p  = mymodel2(fitds,k_m2_dyn+k_m2_dyn_stdv)
Ds_m2_dyn_m  = mymodel2(fitds,k_m2_dyn-k_m2_dyn_stdv)
# Perpendicular
Ds_m2_stat_perp_p = mymodel2(fitds,k_m2_stat_perp+k_m2_stat_perp_stdv)
Ds_m2_stat_perp_m = mymodel2(fitds,k_m2_stat_perp-k_m2_stat_perp_stdv)
Ds_m2_dyn_perp_p  = mymodel2(fitds,k_m2_dyn_perp+k_m2_dyn_perp_stdv)
Ds_m2_dyn_perp_m  = mymodel2(fitds,k_m2_dyn_perp-k_m2_dyn_perp_stdv)

## Model 4
# Parallel
Ds_m4_stat_pk = mymodel4(fitds,k_m4_stat+k_m4_stat_stdv,f_m4_stat)
Ds_m4_stat_mk = mymodel4(fitds,k_m4_stat-k_m4_stat_stdv,f_m4_stat)
Ds_m4_stat_pf = mymodel4(fitds,k_m4_stat,f_m4_stat+f_m4_stat_stdv)
Ds_m4_stat_mf = mymodel4(fitds,k_m4_stat,f_m4_stat-f_m4_stat_stdv)
Ds_m4_dyn_pk  = mymodel4(fitds,k_m4_dyn+k_m4_dyn_stdv,f_m4_dyn)
Ds_m4_dyn_mk  = mymodel4(fitds,k_m4_dyn-k_m4_dyn_stdv,f_m4_dyn)
Ds_m4_dyn_pf  = mymodel4(fitds,k_m4_dyn,f_m4_dyn+f_m4_dyn_stdv)
Ds_m4_dyn_mf  = mymodel4(fitds,k_m4_dyn,f_m4_dyn-f_m4_dyn_stdv)
# Perpendicular
Ds_m4_stat_perp_pk = mymodel4(fitds,k_m4_stat_perp+k_m4_stat_perp_stdv,f_m4_stat_perp)
Ds_m4_stat_perp_mk = mymodel4(fitds,k_m4_stat_perp-k_m4_stat_perp_stdv,f_m4_stat_perp)
Ds_m4_stat_perp_pf = mymodel4(fitds,k_m4_stat_perp,f_m4_stat_perp+f_m4_stat_perp_stdv)
Ds_m4_stat_perp_mf = mymodel4(fitds,k_m4_stat_perp,f_m4_stat_perp-f_m4_stat_perp_stdv)
Ds_m4_dyn_perp_pk  = mymodel4(fitds,k_m4_dyn_perp+k_m4_dyn_perp_stdv,f_m4_dyn_perp)
Ds_m4_dyn_perp_mk  = mymodel4(fitds,k_m4_dyn_perp-k_m4_dyn_perp_stdv,f_m4_dyn_perp)
Ds_m4_dyn_perp_pf  = mymodel4(fitds,k_m4_dyn_perp,f_m4_dyn_perp+f_m4_dyn_perp_stdv)
Ds_m4_dyn_perp_mf  = mymodel4(fitds,k_m4_dyn_perp,f_m4_dyn_perp-f_m4_dyn_perp_stdv)

## Power law 3
# Parallel
Ds_p3_stat_pn = powerlaw3(fitds,n_p3_stat+n_p3_stat_stdv,k_p3_stat)
Ds_p3_stat_mn = powerlaw3(fitds,n_p3_stat-n_p3_stat_stdv,k_p3_stat)
Ds_p3_stat_pk = powerlaw3(fitds,n_p3_stat,k_p3_stat+k_p3_stat_stdv)
Ds_p3_stat_mk = powerlaw3(fitds,n_p3_stat,k_p3_stat-k_p3_stat_stdv)
Ds_p3_dyn_pn  = powerlaw3(fitds,n_p3_dyn+n_p3_dyn_stdv,k_p3_dyn)
Ds_p3_dyn_mn  = powerlaw3(fitds,n_p3_dyn-n_p3_dyn_stdv,k_p3_dyn)
Ds_p3_dyn_pk  = powerlaw3(fitds,n_p3_dyn,k_p3_dyn+k_p3_dyn_stdv)
Ds_p3_dyn_mk  = powerlaw3(fitds,n_p3_dyn,k_p3_dyn-k_p3_dyn_stdv)
# Perpendicular
Ds_p3_stat_perp_pn = powerlaw3(fitds,n_p3_stat_perp+n_p3_stat_perp_stdv,k_p3_stat_perp)
Ds_p3_stat_perp_mn = powerlaw3(fitds,n_p3_stat_perp-n_p3_stat_perp_stdv,k_p3_stat_perp)
Ds_p3_stat_perp_pk = powerlaw3(fitds,n_p3_stat_perp,k_p3_stat_perp+k_p3_stat_perp_stdv)
Ds_p3_stat_perp_mk = powerlaw3(fitds,n_p3_stat_perp,k_p3_stat_perp-k_p3_stat_perp_stdv)
Ds_p3_dyn_perp_pn  = powerlaw3(fitds,n_p3_dyn_perp+n_p3_dyn_perp_stdv,k_p3_dyn_perp)
Ds_p3_dyn_perp_mn  = powerlaw3(fitds,n_p3_dyn_perp-n_p3_dyn_perp_stdv,k_p3_dyn_perp)
Ds_p3_dyn_perp_pk  = powerlaw3(fitds,n_p3_dyn_perp,k_p3_dyn_perp+k_p3_dyn_perp_stdv)
Ds_p3_dyn_perp_mk  = powerlaw3(fitds,n_p3_dyn_perp,k_p3_dyn_perp-k_p3_dyn_perp_stdv)


plt.figure(figsize=(14,5))
ax = plt.subplot(111)
ax.plot(spacings_dyn, Dparallel_dyn, color='g', label=r'$D_\parallel$, dyn.')
ax.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
ax.plot(fitds,Ds_m2_dyn_m, '-.', label=r'$k+\Delta k$') #label=r'k=%.2e' % (k_m2_dyn+k_m2_dyn_stdv))
ax.plot(fitds,Ds_m2_dyn, '--', label=r'$k$')#=%.2e' % k_m2_dyn)
ax.plot(fitds,Ds_m2_dyn_p, '-.', label=r'$k-\Delta k$') #, label=r'k=%.2e' % (k_m2_dyn-k_m2_dyn_stdv))
plt.xlabel(r'$d$ (nm)', fontsize=14)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=14)
plt.title('Custom model 2')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)


plt.figure(figsize=(14,5))
ax = plt.subplot(111)
ax.plot(spacings_dyn, Dparallel_dyn, color='g', label=r'$D_\parallel$, dyn.')
ax.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
ax.plot(fitds,Ds_m4_dyn_mk, '-.', label=r'$k-\Delta k$, $f$') #, label=r'k=%.2e, f=%.2e' % (k_m4_dyn-k_m4_dyn_stdv,f_m4_dyn))
ax.plot(fitds,Ds_m4_dyn_mf, ':', label=r'$k$, $f-\Delta f$') #, label=r'k=%.2e, f=%.2e' % (k_m4_dyn,f_m4_dyn-f_m4_dyn_stdv))
ax.plot(fitds,Ds_m4_dyn, '--', label=r'$k$, $f$') #, label=r'k=%.2e, f=%.2e' % (k_m4_dyn,f_m4_dyn))
ax.plot(fitds,Ds_m4_dyn_pk, '-.', label=r'$k+\Delta k$, $f$')#, label=r'k=%.2e, f=%.2e' % (k_m4_dyn+k_m4_dyn_stdv,f_m4_dyn))
ax.plot(fitds,Ds_m4_dyn_pf, ':', label=r'$k$, $f+\Delta f$')#, label=r'k=%.2e, f=%.2e' % (k_m4_dyn,f_m4_dyn+f_m4_dyn_stdv))
plt.xlabel(r'$d$ (nm)', fontsize=14)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=14)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.title('Custom model 4')


## Plotting:# Need fancy plotting. Find code and redo fancy plotting.
plt.figure(figsize=(14,5))
ax = plt.subplot(111)
ax.plot(spacings_stat, Dparallel_stat, color='limegreen',label=r'$D_\parallel$, stat.')
ax.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen', alpha=0.2)
ax.plot(fitds,Ds_m4_stat_mk, '-.', label=r'$k-\Delta k$, $f$') #, label=r'k=%.2e, f=%.2e' % (k_m4_stat-k_m4_stat_stdv,f_m4_stat))
ax.plot(fitds,Ds_m4_stat_mf, ':', label=r'$k$, $f-\Delta f$') #, label=r'k=%.2e, f=%.2e' % (k_m4_stat,f_m4_stat-f_m4_stat_stdv))
ax.plot(fitds,Ds_m4_stat, '--', label=r'$k$, $f$') #, label=r'k=%.2e, f=%.2e' % (k_m4_stat,f_m4_stat))
ax.plot(fitds,Ds_m4_stat_pk, '-.', label=r'$k+\Delta k$, $f$') #, label=r'k=%.2e, f=%.2e' % (k_m4_stat+k_m4_stat_stdv,f_m4_stat))
ax.plot(fitds,Ds_m4_stat_pf, ':', label=r'$k$, $f+\Delta k$') #, label=r'k=%.2e, f=%.2e' % (k_m4_stat,f_m4_stat+f_m4_stat_stdv))
plt.xlabel(r'$d$ (nm)', fontsize=14)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=14)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.title('Custom model 4')

#... Do the same, men for perp.
plt.figure(figsize=(14,5))
ax = plt.subplot(111)
ax.plot(spacings_dyn, Dzs_dyn, color='g', label=r'$D_\perp$, dyn.')
ax.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='g', alpha=0.2)
ax.plot(fitds,Ds_p3_dyn_perp_mn, '-.', label=r'$n-\Delta n$, $k$') #, label=r'n=%.2e, k=%.2e' % (n_p3_dyn_perp-n_p3_dyn_perp_stdv,k_p3_dyn_perp))
ax.plot(fitds,Ds_p3_dyn_perp_mk, ':', label=r'$n$, $k-\Delta k$') #, label=r'n=%.2e, k=%.2e' % (n_p3_dyn_perp,k_p3_dyn_perp-k_p3_dyn_perp_stdv))
ax.plot(fitds,Ds_p3_dyn_perp, '--', label=r'$n$, $k$') #, label=r'n=%.2e, k=%.2e' % (n_p3_dyn_perp,k_p3_dyn_perp))
ax.plot(fitds,Ds_p3_dyn_perp_pn, '-.', label=r'$n+\Delta n$, $k$') #, label=r'n=%.2e, k=%.2e' % (n_p3_dyn_perp+n_p3_dyn_perp_stdv,k_p3_dyn_perp))
ax.plot(fitds,Ds_p3_dyn_perp_pk, ':', label=r'$n$, $k+\Delta k$') #, label=r'n=%.2e, k=%.2e' % (n_p3_dyn_perp,k_p3_dyn_perp+k_p3_dyn_perp_stdv))
plt.xlabel(r'$d$ (nm)', fontsize=14)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=14)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.title('Power law 3')


plt.figure(figsize=(14,5))
ax = plt.subplot(111)
ax.errorbar(spacings_stat, Dzs_stat, color='limegreen', label=r'$D_\perp$, stat.')
ax.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='limegreen', alpha=0.2)
#ax.plot(fitds,Ds_p3_stat_perp, '-.', label=r'Power law 3')
ax.plot(fitds,Ds_p3_stat_perp_mn, '-.', label=r'$n-\Delta n$, $k$') #, label=r'n=%.2e, k=%.2e' % (n_p3_stat_perp-n_p3_stat_perp_stdv,k_p3_stat_perp))
ax.plot(fitds,Ds_p3_stat_perp_mk, ':', label=r'$n$, $k-\Delta k$') #, label=r'n=%.2e, k=%.2e' % (n_p3_stat_perp,k_p3_stat_perp-k_p3_stat_perp_stdv))
ax.plot(fitds,Ds_p3_stat_perp, '--', label=r'$n$, $k$') #, label=r'n=%.2e, p=%.2e' % (n_p3_stat_perp,p_k3_stat_perp))
ax.plot(fitds,Ds_p3_stat_perp_pn, '-.', label=r'$n+\Delta n$, $k$') #, label=r'n=%.2e, k=%.2e' % (n_p3_stat_perp+n_p3_stat_perp_stdv,k_p3_stat_perp))
ax.plot(fitds,Ds_p3_stat_perp_pk, ':', label=r'$n$, $k+\Delta k$') #, label=r'n=%.2e, k=%.2e' % (n_p3_stat_perp,k_p3_stat_perp+k_p3_stat_perp_stdv))
plt.xlabel(r'$d$ (nm)', fontsize=14)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=14)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.title('Power law 3')
#'''



print('\n\nModel 1')
print('k stat, par:', k_m1_stat, '+/-', k_m1_stat_stdv)
print('k stat, perp:', k_m1_stat_perp, '+/-', k_m1_stat_perp_stdv)
print('k dyn, par:', k_m1_dyn, '+/-', k_m1_dyn_stdv)
print('k dyn, perp:', k_m1_dyn_perp, '+/-', k_m1_dyn_perp_stdv)

print('\n\nModel 2')
print('k stat, par:', k_m2_stat, '+/-', k_m2_stat_stdv)
print('k stat, perp:', k_m2_stat_perp, '+/-', k_m2_stat_perp_stdv)
print('k dyn, par:', k_m2_dyn, '+/-', k_m2_dyn_stdv)
print('k dyn, perp:', k_m2_dyn_perp, '+/-', k_m2_dyn_perp_stdv)


# Model 3
print('\n\nModel 3')
print('k stat, par:', k_m3_stat, '+/-', k_m3_stat_stdv)
print('k stat, perp:', k_m3_stat_perp, '+/-', k_m3_stat_perp_stdv)
print('k dyn, par:', k_m3_dyn, '+/-', k_m3_dyn_stdv)
print('k dyn, perp:', k_m3_dyn_perp, '+/-', k_m3_dyn_perp_stdv)

#'''
# Model 4
print('\n\nModel 4')
print('k stat, par:', k_m4_stat, '+/-', k_m4_stat_stdv,'; f stat, par:', f_m4_stat, '+/-', f_m4_stat_stdv)
print('k stat, perp:', k_m4_stat_perp, '+/-', k_m4_stat_perp_stdv,'; f stat, perp:', f_m4_stat_perp, '+/-', f_m4_stat_perp_stdv)
print('k dyn, par:', k_m4_dyn, '+/-', k_m4_dyn_stdv,'; f dyn, par:', f_m4_dyn, '+/-', f_m4_dyn_stdv)
print('k dyn, perp:', k_m4_dyn_perp, '+/-', k_m4_dyn_perp_stdv,'; f dyn, perp:', f_m4_dyn_perp, '+/-', f_m4_dyn_perp_stdv)
#'''

# Model 6
print('\n\nModel 6')
print('k stat, par:', k_m6_stat, '+/-', k_m6_stat_stdv)
print('k stat, perp:', k_m6_stat_perp, '+/-', k_m6_stat_perp_stdv)
print('k dyn, par:', k_m6_dyn, '+/-', k_m6_dyn_stdv)
print('k dyn, perp:', k_m6_dyn_perp, '+/-', k_m6_dyn_perp_stdv)
print('  ')


# Model 7
print('\n\nModel 7')
print('k stat, par:', k_m7_stat, '+/-', k_m7_stat_stdv)
print('k stat, perp:', k_m7_stat_perp, '+/-', k_m7_stat_perp_stdv)
print('k dyn, par:', k_m7_dyn, '+/-', k_m7_dyn_stdv)
print('k dyn, perp:', k_m7_dyn_perp, '+/-', k_m7_dyn_perp_stdv)
print('  ')

# Model 8
print('\n\nModel 8')
print('k stat, par:', k_m8_stat, '+/-', k_m8_stat_stdv)
print('k stat, perp:', k_m8_stat_perp, '+/-', k_m8_stat_perp_stdv)
print('k dyn, par:', k_m8_dyn, '+/-', k_m8_dyn_stdv)
print('k dyn, perp:', k_m8_dyn_perp, '+/-', k_m8_dyn_perp_stdv)
print('  ')

# Power law 1
print('\n\nPower law 1')
print('k stat, par:', k_p1_stat, '+/-', k_p1_stat_stdv,'l stat, par:', l_p1_stat, '+/-', l_p1_stat_stdv, 'n stat, par:', n_p1_stat, '+/-', n_p1_stat_stdv)
print('k stat, perp:', k_p1_stat_perp, '+/-', k_p1_stat_perp_stdv,'l stat, perp:', l_p1_stat_perp, '+/-', l_p1_stat_perp_stdv, 'n stat, perp:', n_p1_stat_perp, '+/-', n_p1_stat_perp_stdv)
print('k dyn, par:', k_p1_dyn, '+/-', k_p1_dyn_stdv,'l dyn, par:', l_p1_dyn, '+/-', l_p1_dyn_stdv, 'n dyn, par:', n_p1_dyn, '+/-', n_p1_dyn_stdv)
print('k dyn, perp:', k_p1_dyn_perp, '+/-', k_p1_dyn_perp_stdv,'l dyn, perp:', l_p1_dyn_perp, '+/-', l_p1_dyn_perp_stdv, 'n dyn, perp:', n_p1_dyn_perp, '+/-', n_p1_dyn_perp_stdv)
print('  ')

# Power law 2
print('\n\nPower law 2')
print('l stat, par:', l_p2_stat, '+/-', l_p2_stat_stdv, 'n stat, par:', n_p2_stat, '+/-', n_p2_stat_stdv)
print('l stat, perp:', l_p2_stat_perp, '+/-', l_p2_stat_perp_stdv, 'n stat, perp:', n_p2_stat_perp, '+/-', n_p2_stat_perp_stdv)
print('l dyn, par:', l_p2_dyn, '+/-', l_p2_dyn_stdv, 'n dyn, par:', n_p2_dyn, '+/-', n_p2_dyn_stdv)
print('l dyn, perp:', l_p2_dyn_perp, '+/-', l_p2_dyn_perp_stdv, 'n dyn, perp:', n_p2_dyn_perp, '+/-', n_p2_dyn_perp_stdv)
print('  ')


# Power law 3
print('\n\nPower law 3')
print('n stat, par:', n_p3_stat, '+/-', n_p3_stat_stdv, 'k stat, par:', k_p3_stat, '+/-', k_p3_stat_stdv)
print('n stat, perp:', n_p3_stat_perp, '+/-', n_p3_stat_perp_stdv, 'k stat, perp:', k_p3_stat_perp, '+/-', n_p3_stat_perp_stdv)
print('n dyn, par:', n_p3_dyn, '+/-', n_p3_dyn_stdv, 'k dyn, par:', k_p3_dyn, '+/-', k_p3_dyn_stdv)
print('n dyn, perp:', n_p3_dyn_perp, '+/-', n_p3_dyn_perp_stdv, 'k dyn, perp:', k_p3_dyn_perp, '+/-', n_p3_dyn_perp_stdv)
print('  ')

# Power law 4
print('\n\nPower law 4')
print('n stat, par:', n_p4_stat, '+/-', n_p4_stat_stdv)
print('n stat, perp:', n_p4_stat_perp, '+/-', n_p4_stat_perp_stdv)
print('n dyn, par:', n_p4_dyn, '+/-', n_p4_dyn_stdv)
print('n dyn, perp:', n_p4_dyn_perp, '+/-', n_p4_dyn_perp_stdv)
print('  ')


### From Shen & Chen:

print('-----------------------')
print('Ordered packings:')
print('sigma, op, dyn, par:',sigma_op_dyn,'+/-',sigma_op_dyn_rms)
print('sigma, op, stat, par:',sigma_op_stat,'+/-',sigma_op_stat_rms)
print('sigma, op, dyn, perp:',sigma_op_dyn_perp,'+/-',sigma_op_dyn_perp_rms)
print('sigma, op, stat, perp:',sigma_op_stat_perp,'+/-',sigma_op_stat_perp_rms)
print('-----------------------')
print('Hyperbola of revolution:')
print('sigma, hr, dyn, par:',sigma_hr_dyn,'+/-',sigma_hr_dyn_rms)
print('sigma, hr, stat, par:',sigma_hr_stat,'+/-',sigma_hr_stat_rms)
print('sigma, hr, dyn, perp:',sigma_hr_dyn_perp,'+/-',sigma_hr_dyn_perp_rms)
print('sigma, hr, stat, perp:',sigma_hr_stat_perp,'+/-',sigma_hr_stat_perp_rms)
print('-----------------------')
print('Heterogeneous catalyst:')
print('sigma, hc, dyn, par:',sigma_hc_dyn,'+/-',sigma_hc_dyn_rms)
print('sigma, hc, stat, par:',sigma_hc_stat,'+/-',sigma_hc_stat_rms)
print('sigma, hc, dyn, perp:',sigma_hc_dyn_perp,'+/-',sigma_hc_dyn_perp_rms)
print('sigma, hc, stat, perp:',sigma_hc_stat_perp,'+/-',sigma_hc_stat_perp_rms)
print('-----------------------')
print('Cation-exchange resin membrane:')
print('sigma, rm, dyn, par:',sigma_rm_dyn,'+/-',sigma_rm_dyn_rms)
print('sigma, rm, stat, par:',sigma_rm_stat,'+/-',sigma_rm_stat_rms)
print('sigma, rm, dyn, perp:',sigma_rm_dyn_perp,'+/-',sigma_rm_dyn_perp_rms)
print('sigma, rm, stat, perp:',sigma_rm_stat_perp,'+/-',sigma_rm_stat_perp_rms)
print('-----------------------')



''' ## Made some changes, can't do this anymore
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
# Cation exchange resin membrane 8
D_rmsd_rm_rw   = rmsd(DRW,Ds_rm_rw)
D_rmsd_rm_stat = rmsd(Dparallel_stat,Ds_rm_stat)
D_rmsd_rm_dyn  = rmsd(Dparallel_dyn,Ds_rm_dyn)
D_rmsd_rm_f    = rmsd(Dparallel_f,Ds_rm_f)


# Write rmsd to file

#rmsdfilename_RW
#rmsdfilename_forest
#rmsdfilename_dyn
#rmsdfilename_stat


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

#####
rmsdfile_rw.write('%.16f %.16f %.16f %.16f\n' % (D_rmsd_op_rw,D_rmsd_hr_rw,D_rmsd_hc_rw,D_rmsd_rm_rw))
rmsdfile_rw.write('sigma: %.16f %.16f %.16f %.16f\n' % (sigma_op_rw,sigma_hr_rw,sigma_hc_rw,sigma_rm_rw))
rmsdfile_rw.close()

rmsdfile_f = open(rmsdfilename_forest,'w')
rmsdfile_f.write('%.16f %.16f %.16f %.16f\n' % (D_rmsd_op_f,D_rmsd_hr_f,D_rmsd_hc_f,D_rmsd_rm_f))
rmsdfile_f.write('sigma: %.16f %.16f %.16f %.16f\n' % (sigma_op_f,sigma_hr_f,sigma_hc_f,sigma_rm_f))
rmsdfile_f.close()

rmsdfile_dyn = open(rmsdfilename_dyn,'w')
rmsdfile_dyn.write('%.16f %.16f %.16f %.16f\n' % (D_rmsd_op_dyn,D_rmsd_hr_dyn,D_rmsd_hc_dyn,D_rmsd_rm_dyn))
rmsdfile_dyn.write('sigma: %.16f %.16f %.16f %.16f\n' % (sigma_op_dyn,sigma_hr_dyn,sigma_hc_dyn,sigma_rm_dyn))
rmsdfile_dyn.close()

rmsdfile_stat = open(rmsdfilename_stat,'w')
rmsdfile_stat.write('%.16f %.16f %.16f %.16f\n' % (D_rmsd_op_stat,D_rmsd_hr_stat,D_rmsd_hc_stat,D_rmsd_rm_stat))
rmsdfile_stat.write('sigma: %.16f %.16f %.16f %.16f\n' % (sigma_op_stat,sigma_hr_stat,sigma_hc_stat,sigma_rm_stat))
rmsdfile_stat.close()
'''

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

print('Covariance matrices:')
print('Modell 4:')
print('pcov4_statpar:',pcov4_statpar)
print('pcov4_dynpar:',pcov4_dynpar)
print('pcov4_statperp:',pcov4_statperp)
print('pcov4_dynperp:',pcov4_dynperp)

print('Power law 3:')
print('pcovp3_statpar:',pcovp3_statpar)
print('pcovp3_dynpar:',pcovp3_dynpar)
print('pcovp3_statperp:',pcovp3_statperp)
print('pcovp3_dynperp:',pcovp3_dynperp)

'''
plotname_sigma_bare_dynamic          = endlocation + 'D_vs_d_varioussigmas_dynamic.png'
plotname_sigma_bare_static           = endlocation + 'D_vs_d_varioussigmas_static.png' 
plotname_difftime_sigma_bare_dynamic = endlocation + 'diffusiontime_dynamic.png'
plotname_difftime_sigma_bare_static  = endlocation + 'diffusiontime_static.png' 
# Loop and print graphs?
thickness = 50e-9 # m
sigmas = np.array([0.1,0.5,1.0,1.5,2.0,3.0])
Nsigmas = len(sigmas)
fitds_dyn_store = []
diffusiontimes_dyn = []
plt.figure(figsize=(6,5))
for sigma in sigmas:
    Ds_p3_dyn_perp = powerlaw3_sigmas(fitds,sigma,n_p3_dyn_perp,k_p3_dyn_perp)
    fitds_II          = []
    Ds_p3_dyn_perp_II = []
    diffusiontime_dyn = []
    for i in range(len(Ds_p3_dyn_perp)):
        if Ds_p3_dyn_perp[i]>0:
            thisD = Ds_p3_dyn_perp[i]
            fitds_II.append(fitds[i])
            Ds_p3_dyn_perp_II.append(thisD)
            D = thisD*Dzs_bulk
            difftime = thickness**2/(2.*D)
            diffusiontime_dyn.append(difftime)
    plt.plot(fitds_II,Ds_p3_dyn_perp_II,label=r'$\sigma_{\mathregular{LJ}}$ = %s nm' %str(sigma))
    diffusiontimes_dyn.append(diffusiontime_dyn)
    fitds_dyn_store.append(fitds_II)
plt.xlabel('d (nm)')
plt.ylabel(r'$D_\perp/D_{\mathregular{bulk}}$, dynamic brush')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(plotname_sigma_bare_dynamic)

fitds_stat_store = []
diffusiontimes_stat = []
plt.figure(figsize=(6,5))
for sigma in sigmas:
    Ds_p3_stat_perp = powerlaw3_sigmas(fitds,sigma,n_p3_stat_perp,k_p3_stat_perp)
    fitds_II           = []
    Ds_p3_stat_perp_II = []
    diffusiontime_stat = []
    for i in range(len(Ds_p3_stat_perp)):
        if Ds_p3_stat_perp[i]>0:
            thisD = Ds_p3_stat_perp[i]
            fitds_II.append(fitds[i])
            Ds_p3_stat_perp_II.append(thisD)
            D = thisD*Dzs_bulk
            difftime = thickness**2/(2.*D)
            diffusiontime_stat.append(difftime)
    plt.plot(fitds_II,Ds_p3_stat_perp_II,label=r'$\sigma_{\mathregular{LJ}}$ = %s nm' %str(sigma))
    diffusiontimes_stat.append(diffusiontime_stat)
    fitds_stat_store.append(fitds_II)
plt.xlabel('d (nm)')
plt.ylabel(r'$D_\perp/D_{\mathregular{bulk}}$, static brush')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(plotname_sigma_bare_static)


plt.figure(figsize=(6,5))
for i in range(Nsigmas):
    plt.plot(fitds_dyn_store[i],diffusiontimes_dyn[i],label=r'$\sigma_{\mathregular{LJ}}$ = %s nm' %str(sigmas[i]))
plt.xlabel('d (nm)')
plt.ylabel('Diffusion time (s)')
plt.title('Net thickness %s nm, dynamic brush' % str(thickness*1e9)) # Or change to nm?
plt.legend(loc='upper right')
plt.savefig(plotname_difftime_sigma_bare_dynamic)

plt.figure(figsize=(6,5))
for i in range(Nsigmas):
    plt.plot(fitds_stat_store[i],diffusiontimes_stat[i],label=r'$\sigma_{\mathregular{LJ}}$ = %s nm' %str(sigmas[i]))
plt.xlabel('d (nm)')
plt.ylabel('Diffusion time (s)')
plt.title('Net thickness %s nm, static brush' % str(thickness*1e9)) # Or change to nm?
plt.legend(loc='upper right')
plt.savefig(plotname_difftime_sigma_bare_static)

### sigma_walker: #############################################################################################################################
plotname_sigma_free_dynamic          = endlocation + 'D_vs_d_varioussigmas_sigmafree_dynamic.png'
plotname_sigma_free_static           = endlocation + 'D_vs_d_varioussigmas_sigmafree_static.png' 
plotname_difftime_sigma_free_dynamic = endlocation + 'diffusiontime_dynamic_sigmafree.png'
plotname_difftime_sigma_free_static  = endlocation + 'diffusiontime_static_sigmafree.png' 

thickness = 50e-9 # m
sigma_walker = np.array([0.01,0.1,0.5,1.0,1.5,2.0,3.0])
sigma_chain  = 1.0
sigmas = (sigma_walker+sigma_chain)/2.#np.array([0.1,0.5,1.0,1.5,2.0,3.0])
Nsigmas = len(sigmas)
fitds_dyn_store = []
diffusiontimes_dyn = []
plt.figure(figsize=(6,5))
for sigma in sigma_walker:
    Ds_p3_dyn_perp = powerlaw3_sigmas(fitds,sigma,n_p3_dyn_perp,k_p3_dyn_perp)
    fitds_II          = []
    Ds_p3_dyn_perp_II = []
    diffusiontime_dyn = []
    for i in range(len(Ds_p3_dyn_perp)):
        if Ds_p3_dyn_perp[i]>0:
            thisD = Ds_p3_dyn_perp[i]
            fitds_II.append(fitds[i])
            Ds_p3_dyn_perp_II.append(thisD)
            D = thisD*Dzs_bulk
            difftime = thickness**2/(2.*D)
            diffusiontime_dyn.append(difftime)
    plt.plot(fitds_II,Ds_p3_dyn_perp_II,label=r'$\sigma_{\mathregular{LJ}}$ = %s nm' %str(sigma))
    diffusiontimes_dyn.append(diffusiontime_dyn)
    fitds_dyn_store.append(fitds_II)
plt.xlabel('d (nm)')
plt.ylabel(r'$D_\perp/D_{\mathregular{bulk}}$, dynamic brush')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(plotname_sigma_free_dynamic)

fitds_stat_store = []
diffusiontimes_stat = []
plt.figure(figsize=(6,5))
for sigma in sigma_walker:
    Ds_p3_stat_perp = powerlaw3_sigmas(fitds,sigma,n_p3_stat_perp,k_p3_stat_perp)
    fitds_II           = []
    Ds_p3_stat_perp_II = []
    diffusiontime_stat = []
    for i in range(len(Ds_p3_stat_perp)):
        if Ds_p3_stat_perp[i]>0:
            thisD = Ds_p3_stat_perp[i]
            fitds_II.append(fitds[i])
            Ds_p3_stat_perp_II.append(thisD)
            D = thisD*Dzs_bulk
            difftime = thickness**2/(2.*D)
            diffusiontime_stat.append(difftime)
    plt.plot(fitds_II,Ds_p3_stat_perp_II,label=r'$\sigma_{\mathregular{LJ}}$ = %s nm' %str(sigma))
    diffusiontimes_stat.append(diffusiontime_stat)
    fitds_stat_store.append(fitds_II)
plt.xlabel('d (nm)')
plt.ylabel(r'$D_\perp/D_{\mathregular{bulk}}$, static brush')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(plotname_sigma_free_static)


plt.figure(figsize=(6,5))
for i in range(Nsigmas):
    plt.plot(fitds_dyn_store[i],diffusiontimes_dyn[i],label=r'$\sigma_{\mathregular{LJ}}$ = %s nm' %str(sigma_walker[i]))
plt.xlabel('d (nm)')
plt.ylabel('Diffusion time (s)')
plt.title('Net thickness %s nm, dynamic brush' % str(thickness*1e9)) # Or change to nm?
plt.legend(loc='upper right')
plt.savefig(plotname_difftime_sigma_free_dynamic)

plt.figure(figsize=(6,5))
for i in range(Nsigmas):
    plt.plot(fitds_stat_store[i],diffusiontimes_stat[i],label=r'$\sigma_{\mathregular{LJ}}$ = %s nm' %str(sigma_walker[i]))
plt.xlabel('d (nm)')
plt.ylabel('Diffusion time (s)')
plt.title('Net thickness %s nm, static brush' % str(thickness*1e9)) # Or change to nm?
plt.legend(loc='upper right')
plt.savefig(plotname_difftime_sigma_free_static)
'''

### sigma_walker, scaling the mass : #############################################################################################################################
plotname_sigma_free_massadj_dynamic                  = endlocation + 'D_vs_d_varioussigmas_sigmafree_dynamic_adjustingmass.png'
plotname_sigma_free_massadj_static                   = endlocation + 'D_vs_d_varioussigmas_sigmafree_static_adjustingmass.png' 
plotname_difftime_sigma_free_massadj_dynamic         = endlocation + 'diffusiontime_dynamic_sigmafree_adjustingmass.png'
plotname_difftime_sigma_free_massadj_static          = endlocation + 'diffusiontime_static_sigmafree_adjustingmass.png' 
plotname_difftime_sigma_free_massadj_dynamic_zooming = endlocation + 'diffusiontime_dynamic_sigmafree_adjustingmass_zoomzoom.png'
plotname_difftime_sigma_free_massadj_static_zooming  = endlocation + 'diffusiontime_static_sigmafree_adjustingmass_zoomzoom.png' 

bulkfilename_II = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Bulk/Varysigmas/D_vs_sigma_better_rms_Nestimates10.txt'
bulkfile_II = open(bulkfilename_II,'r')
hline = bulkfile_II.readline()
lines = bulkfile_II.readlines()

sigma_walker = []
Dzs_bulk_vm = []
for line in lines:
    words = line.split()
    if len(words)>0:
        sigma_walker.append(float(words[0]))
        Dzs_bulk_vm.append(float(words[3]))
bulkfile_II.close()

mylinestyles = ['-','--','-.','-,',':','-']

print('sigma_walker:',sigma_walker)

thickness = 50e-9 # m
sigma_chain  = 1.0
sigma_walker = np.array(sigma_walker) # np.array([0.1,0.5,1.0,1.5,2.0,3.0,5.0,10.0])
sigmas = (sigma_walker+sigma_chain)/2.
Nsigmas = len(sigmas)
fitds_dyn_store = []
diffusiontimes_dyn = []
j = 0
plt.figure(figsize=(7.5,5))
#ax = plt.subplot(111)
for sigma in sigma_walker:
    Ds_p3_dyn_perp = powerlaw3_sigmas(fitds,sigma,n_p3_dyn_perp,k_p3_dyn_perp)
    fitds_II          = []
    Ds_p3_dyn_perp_II = []
    diffusiontime_dyn = []
    for i in range(len(Ds_p3_dyn_perp)):
        if Ds_p3_dyn_perp[i]>0:
            thisD = Ds_p3_dyn_perp[i]
            fitds_II.append(fitds[i])
            Ds_p3_dyn_perp_II.append(thisD)
            D = thisD*Dzs_bulk_vm[j]
            difftime = thickness**2/(2.*D)
            diffusiontime_dyn.append(difftime)
    plt.plot(fitds_II,Ds_p3_dyn_perp_II,mylinestyles[j],label=r'$\sigma_{\mathregular{LJ}}$ = %s nm' %str(sigma))
    diffusiontimes_dyn.append(diffusiontime_dyn)
    fitds_dyn_store.append(fitds_II)
    j+=1
plt.xlabel('d (nm)', fontsize=14)
plt.ylabel(r'$D_\perp/D_{\mathregular{bulk}}$', fontsize=14)#, dynamic brush', fontsize=14)
plt.legend(bbox_to_anchor=(1,0), loc="lower left", fontsize=11)
plt.tight_layout()
#plt.legend(loc='lower right')
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_sigma_free_massadj_dynamic)

j = 0 
fitds_stat_store = []
diffusiontimes_stat = []
plt.figure(figsize=(7.5,5))
#ax = plt.subplot(111)
for sigma in sigma_walker:
    Ds_p3_stat_perp = powerlaw3_sigmas(fitds,sigma,n_p3_stat_perp,k_p3_stat_perp)
    fitds_II           = []
    Ds_p3_stat_perp_II = []
    diffusiontime_stat = []
    for i in range(len(Ds_p3_stat_perp)):
        if Ds_p3_stat_perp[i]>0:
            thisD = Ds_p3_stat_perp[i]
            fitds_II.append(fitds[i])
            Ds_p3_stat_perp_II.append(thisD)
            D = thisD*Dzs_bulk_vm[j]
            difftime = thickness**2/(2.*D)
            diffusiontime_stat.append(difftime)
    plt.plot(fitds_II,Ds_p3_stat_perp_II,mylinestyles[j],label=r'$\sigma_{\mathregular{LJ}}$ = %s nm' %str(sigma))
    diffusiontimes_stat.append(diffusiontime_stat)
    fitds_stat_store.append(fitds_II)
    j+=1
plt.xlabel('d (nm)', fontsize=14)
plt.ylabel(r'$D_\perp/D_{\mathregular{bulk}}$', fontsize=14)#, static brush', fontsize=14)
plt.legend(bbox_to_anchor=(1,0), loc="lower left", fontsize=11)
plt.tight_layout()
#plt.legend(loc='lower right')
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_sigma_free_massadj_static)

'''
plt.figure(figsize=(6,5))
for i in range(Nsigmas):
    plt.plot(fitds_dyn_store[i],diffusiontimes_dyn[i],label=r'$\sigma_{\mathregular{LJ}}$ = %s nm' %str(sigma_walker[i]))
plt.xlabel('d (nm)', fontsize=14)
plt.ylabel('Diffusion time (s)', fontsize=14)
#plt.title('Net thickness %s nm, dynamic brush, adjusting mass' % str(thickness*1e9)) # Or change to nm?
plt.legend(loc='upper left', fontsize=11)
plt.tight_layout()
plt.savefig(plotname_difftime_sigma_free_massadj_dynamic)

plt.figure(figsize=(6,5))
for i in range(Nsigmas):
    plt.plot(fitds_stat_store[i],diffusiontimes_stat[i],label=r'$\sigma_{\mathregular{LJ}}$ = %s nm' %str(sigma_walker[i]))
plt.xlabel('d (nm)', fontsize=14)
plt.ylabel('Diffusion time (s)', fontsize=14)
#plt.title('Net thickness %s nm, static brush, adjusting mass' % str(thickness*1e9)) # Or change to nm?
plt.legend(loc='upper left', fontsize=11)
plt.tight_layout()
plt.savefig(plotname_difftime_sigma_free_massadj_static)
'''

###### Zoom-friendly: #######
#'''
plt.figure(figsize=(7.5,5))
for i in range(Nsigmas):#(5):#
    plt.plot(fitds_dyn_store[i],diffusiontimes_dyn[i],mylinestyles[i],label=r'$\sigma_{\mathregular{LJ}}$ = %s nm' %str(sigma_walker[i]))
plt.xlabel('d (nm)', fontsize=14)
plt.ylabel('Diffusion time (s)', fontsize=14)
#plt.title('Net thickness %s nm, dynamic brush, adjusting mass' % str(thickness*1e9)) # Or change to nm?
plt.legend(bbox_to_anchor=(1,0), loc="lower left", fontsize=11)
plt.tight_layout()
plt.savefig(plotname_difftime_sigma_free_massadj_dynamic_zooming)
'''

plt.figure(figsize=(7.5,5))
for i in range(Nsigmas):#(5):#
    plt.plot(fitds_stat_store[i],diffusiontimes_stat[i],label=r'$\sigma_{\mathregular{LJ}}$ = %s nm' %str(sigma_walker[i]))
plt.xlabel('d (nm)')
plt.ylabel('Diffusion time (s)')
#plt.title('Net thickness %s nm, static brush, adjusting mass' % str(thickness*1e9)) # Or change to nm?
plt.legend(bbox_to_anchor=(1,0), loc="lower left")
plt.tight_layout()
plt.savefig(plotname_difftime_sigma_free_massadj_static_zooming)
#'''
##################################################################################################################################################################


if showplots==True:
    plt.show()
