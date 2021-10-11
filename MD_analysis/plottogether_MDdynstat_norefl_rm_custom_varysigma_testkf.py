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

def catex_resin_membr(d,sigma):
    phi = 1 - np.pi*(sigma/d)
    tau = ((2.-phi)/phi)**2
    return 1./tau

def mymodel4(d,k,f,sigma):
    # Assumes probability of neighbor being occupied when you're occupied: kvp**f
    phi = 1 - np.pi*(sigma/d)
    vp = 1-phi
    nevner = 1-k*vp**f
    tau = (1+vp/nevner)**2
    return 1./tau


def get_porosities_probnext(d,k,f,sigma):
    # Just a verification that everything makes sense
    phi = 1 - np.pi*(sigma/d)
    vp = 1-phi
    probnext = k*vp**f
    return phi, probnext


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

# Cation exchange resin membrane
popt, pcov           = curve_fit(catex_resin_membr, spacings_stat, Dparallel_stat, maxfev=thismaxfev)
sigma_rm_stat        = popt[0]
sigma_rm_stat_rms    = np.sqrt(np.diag(pcov))[0]
popt, pcov           = curve_fit(catex_resin_membr, spacings_dyn, Dparallel_dyn, maxfev=thismaxfev)
sigma_rm_dyn         = popt[0]
sigma_rm_dyn_rms     = np.sqrt(np.diag(pcov))[0]
#Perp
popt, pcov             = curve_fit(catex_resin_membr, spacings_stat, Dzs_stat, maxfev=thismaxfev)
sigma_rm_stat_perp     = popt[0]
sigma_rm_stat_perp_rms = np.sqrt(np.diag(pcov))[0]
popt, pcov             = curve_fit(catex_resin_membr, spacings_dyn, Dzs_dyn, maxfev=thismaxfev)
sigma_rm_dyn_perp      = popt[0]
sigma_rm_dyn_perp_rms = np.sqrt(np.diag(pcov))[0]

# Ds
Ds_rm_stat = catex_resin_membr(fitds,sigma_rm_stat)
Ds_rm_dyn  = catex_resin_membr(fitds,sigma_rm_dyn)
Ds_rm_stat_perp = catex_resin_membr(fitds,sigma_rm_stat_perp)
Ds_rm_dyn_perp  = catex_resin_membr(fitds,sigma_rm_dyn_perp)

# My models:
### mymodel4(d,k,f)
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
# Ds
Ds_m4_stat = mymodel4(fitds,k_m4_stat,f_m4_stat,sigma_m4_stat)
Ds_m4_dyn  = mymodel4(fitds,k_m4_dyn,f_m4_dyn,sigma_m4_dyn)
Ds_m4_stat_perp = mymodel4(fitds,k_m4_stat_perp,f_m4_stat_perp,sigma_m4_stat_perp)
Ds_m4_dyn_perp  = mymodel4(fitds,k_m4_dyn_perp,f_m4_dyn_perp,sigma_m4_dyn_perp)

#### Plotting:
##params = {'mathtext.default': 'regular' }  # Ditch this.   
###plt.rcParams.update(params)

# Changing params. Need some bool for f eventually
testk = 0.17  ## Change this to test

# Testing k's:
Ds_m4_dyn_testfit = mymodel4(fitds,testk,f_m4_dyn,sigma_m4_dyn)
Ds_m4_dyn_perp_testfit = mymodel4(fitds,testk,f_m4_dyn_perp,sigma_m4_dyn_perp)
Ds_m4_stat_testfit = mymodel4(fitds,testk,f_m4_stat,sigma_m4_stat)
Ds_m4_stat_perp_testfit = mymodel4(fitds,testk,f_m4_stat_perp,sigma_m4_stat_perp)


plt.figure(figsize=(6,5))#(14,5))
ax = plt.subplot(111)
ax.plot(spacings_dyn, Dparallel_dyn, color='g', label=r'$D_\parallel$, dyn.')
ax.fill_between(spacings_dyn, Dparallel_dyn+Dparallel_stdv_dyn, Dparallel_dyn-Dparallel_stdv_dyn, facecolor='g', alpha=0.2)
#line1, = ax.plot(fitds,Ds_rm_dyn, '--,', label=r'Cation-exchange resin membrane') 
#line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
ax.plot(fitds,Ds_m4_dyn, '-.', color='rebeccapurple', label=r'Custom model, k=%.2f' % k_m4_dyn)
ax.plot(fitds,Ds_m4_dyn_testfit,'--', label=r'Custom model, k=%.2f' % testk)
# Add more plotting! # Need to beware of negative porosity for some!
plt.xlabel(r'$d$ (nm)', fontsize=12)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=12)
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.title(r'$D_\parallel$, dyn.')
plt.legend(loc='lower right')
plt.tight_layout()


## Plotting:
plt.figure(figsize=(6,5))#(14,5))
ax = plt.subplot(111)
ax.plot(spacings_stat, Dparallel_stat, color='limegreen',label=r'$D_\parallel$, stat.')
ax.fill_between(spacings_stat, Dparallel_stat+Dparallel_stdv_stat, Dparallel_stat-Dparallel_stdv_stat, facecolor='limegreen', alpha=0.2)
#line1, = ax.plot(fitds,Ds_rm_stat, '--,', label=r'Cation-exchange resin membrane') 
#line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
ax.plot(fitds,Ds_m4_stat, '-.', color='rebeccapurple', label=r'Custom model, k=%.2f' % k_m4_stat)
ax.plot(fitds,Ds_m4_stat_testfit, '--', label=r'Custom model, k=%.2f' % testk)
plt.xlabel(r'$d$ (nm)')
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=12)
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.title(r'$D_\parallel$, stat.')
plt.legend(loc='lower right')
plt.tight_layout()


#... Do the same, men for perp.
plt.figure(figsize=(6,5))#(14,5))
ax = plt.subplot(111)
ax.plot(spacings_dyn, Dzs_dyn, color='g', label=r'$D_\perp$, dyn.')
ax.fill_between(spacings_dyn, Dzs_dyn+Dzs_stdv_dyn, Dzs_dyn-Dzs_stdv_dyn, facecolor='g', alpha=0.2)
#line1, = ax.plot(fitds,Ds_rm_dyn_perp, '--,', label=r'Cation-exchange resin membrane') 
#line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
ax.plot(fitds,Ds_m4_dyn_perp, '-.', color='rebeccapurple', label=r'Custom model, k=%.2f' % k_m4_dyn_perp)
ax.plot(fitds,Ds_m4_dyn_perp_testfit, '--', label=r'Custom model, k=%.2f' % testk)
# Add more plotting! # Need to beware of negative porosity for some!
plt.xlabel(r'$d$ (nm)', fontsize=12)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=12)
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.title(r'$D_\perp$, dyn.')
plt.legend(loc='lower right')
plt.tight_layout()


plt.figure(figsize=(6,5))#(14,5))
ax = plt.subplot(111)
ax.plot(spacings_stat, Dzs_stat, color='limegreen', label=r'$D_\perp$, stat.')
ax.fill_between(spacings_stat, Dzs_stat+Dzs_stdv_stat, Dzs_stat-Dzs_stdv_stat, facecolor='limegreen', alpha=0.2)
#line1, = ax.plot(fitds,Ds_rm_stat_perp, '--,', label=r'Cation-exchange resin membrane') 
#line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
ax.plot(fitds,Ds_m4_stat_perp, '-.', color='rebeccapurple', label=r'Custom model, k=%.2f' % k_m4_stat_perp)
ax.plot(fitds,Ds_m4_stat_perp_testfit, '--', label=r'Custom model, k=%.2f' % testk)
# Add more plotting! # Need to beware of negative porosity for some!
plt.xlabel(r'$d$ (nm)', fontsize=12)
plt.ylabel(r'$D/D_{\mathregular{bulk}}$', fontsize=12)
#plt.title(r'$D/D_{\mathregular{bulk}}$ vs $d$')
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.title(r'$D_\perp$, stat.')
plt.legend(loc='lower right')
plt.tight_layout()


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




# Shen & Chen-models:


print('Cation-exchange resin membrane:')
print('Sigma, stat, par:', sigma_rm_stat, '+/-', sigma_rm_stat_rms)
print('Sigma, stat, perp:', sigma_rm_stat_perp, '+/-', sigma_rm_stat_perp_rms)
print('Sigma, dyn, par:', sigma_rm_dyn, '+/-', sigma_rm_dyn_rms)
print('Sigma, dyn, perp:', sigma_rm_dyn_perp, '+/-', sigma_rm_dyn_perp_rms)
print('   ')

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

print('\n--------------------')
print('testk:',testk)
rmsd_stat_perp = rmsd(Ds_m4_stat_perp,Ds_m4_stat_perp_testfit)
rmsd_stat_par  = rmsd(Ds_m4_stat,Ds_m4_stat_testfit)
rmsd_dyn_perp  = rmsd(Ds_m4_dyn_perp,Ds_m4_dyn_perp_testfit)
rmsd_dyn_par   = rmsd(Ds_m4_dyn,Ds_m4_dyn_testfit)
print('RMSD:')
print('Stat, perp:', rmsd_stat_perp)
print('Stat, par:', rmsd_stat_par)
print('Dyn, perp:', rmsd_dyn_perp)
print('Dyn, par:', rmsd_dyn_par)

if showplots==True:
    plt.show()

testk_par = np.array([0.5,0.6,0.7,0.8,0.85,0.9,0.95,0.96,0.97,0.98,0.99,0.999,1.0,1.001,1.01,1.1])
testk_dynperp = np.array([0.01,0.02,0.03,0.04,0.05,0.06,0.065,0.067,0.068,0.069,0.07,0.071,0.072,0.073,0.074,0.075,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17])
testk_statperp = np.array([0.5,0.6,0.7,0.77,0.8,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,0.97,1.0])

testk_dynpar = np.linspace(0,1.6,1001)
testk_statpar = np.linspace(0,1.2,1001)
testk_dynperp = np.linspace(0,0.7,301)
testk_statperp = np.linspace(0,1,1001)


Npar_dyn   = len(testk_dynpar)
Npar_stat  = len(testk_statpar)
Nperp_dyn  = len(testk_dynperp)
Nperp_stat = len(testk_statperp)

rms_par_dyn   = np.zeros(Npar_dyn)
rms_par_stat  = np.zeros(Npar_stat)
rms_perp_dyn  = np.zeros(Nperp_dyn)
rms_perp_stat = np.zeros(Nperp_stat)


for i in range(Npar_dyn):
    Ds_m4_dyn_testfit  = mymodel4(fitds,testk_dynpar[i],f_m4_dyn,sigma_m4_dyn)
    rms_par_dyn[i]  = rmsd(Ds_m4_dyn,Ds_m4_dyn_testfit)

for i in range(Npar_dyn):
    Ds_m4_stat_testfit = mymodel4(fitds,testk_statpar[i],f_m4_stat,sigma_m4_stat)
    rms_par_stat[i] = rmsd(Ds_m4_stat,Ds_m4_stat_testfit)

for i in range(Nperp_dyn):
    Ds_m4_dyn_perp_testfit = mymodel4(fitds,testk_dynperp[i],f_m4_dyn_perp,sigma_m4_dyn_perp)
    rms_perp_dyn[i]  = rmsd(Ds_m4_dyn_perp,Ds_m4_dyn_perp_testfit)

for i in range(Nperp_stat):
    Ds_m4_stat_perp_testfit = mymodel4(fitds,testk_statperp[i],f_m4_stat_perp,sigma_m4_stat_perp)
    rms_perp_stat[i]  = rmsd(Ds_m4_stat_perp,Ds_m4_stat_perp_testfit)

plt.figure(figsize=(6,5))
plt.plot(testk_dynpar,rms_par_dyn,label=r'$\parallel$, dyn.')
plt.plot(testk_dynperp,rms_perp_dyn,label=r'$\perp$, dyn.')
plt.plot(testk_statpar,rms_par_stat,label=r'$\parallel$, stat.')
plt.plot(testk_statperp,rms_perp_stat,label=r'$\perp$, stat.')
plt.xlabel('k')
plt.ylabel('RMSD between fit and optimum fit')
plt.legend(loc='upper left')
plt.show()

### Verification that everything is ok:
# Dyn
phi_dynpar, probnext_dynpar = get_porosities_probnext(spacings_dyn,k_m4_dyn,f_m4_dyn,sigma_m4_dyn)
phi_dynperp, probnext_dynperp = get_porosities_probnext(spacings_dyn,k_m4_dyn_perp,f_m4_dyn_perp,sigma_m4_dyn_perp)
# Stat
phi_statpar, probnext_statpar = get_porosities_probnext(spacings_stat,k_m4_stat,f_m4_stat,sigma_m4_stat)
phi_statperp, probnext_statperp = get_porosities_probnext(spacings_stat,k_m4_stat_perp,f_m4_stat_perp,sigma_m4_stat_perp)

print('Verifying models:')
print('Dyn, par:')
print('Porosities:', phi_dynpar)
print('Probnext:', probnext_dynpar)
print('------------')
print('Dyn, perp:')
print('Porosities:', phi_dynperp)
print('Probnext:', probnext_dynperp)
print('------------')
print('Stat, par:')
print('Porosities:', phi_statpar)
print('Probnext:', probnext_statpar)
print('------------')
print('Stat, perp:')
print('Porosities:', phi_statperp)
print('Probnext:', probnext_statperp)
print('------------')


ktest = 5.523210923300265e-07
ftest = 0.001078013335684682
sigmatest = 0.9086087538759112
print('k=',ktest)
print('f=',ftest)
print('sigma=',sigmatest)
phi_test, probnext_test = get_porosities_probnext(spacings_stat,ktest,ftest,sigmatest)
print('Porosities;',phi_test)
print('Probnext;',probnext_test)
