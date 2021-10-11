import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
import numpy as np
import random
import math
import time
import os
import glob

def myexponential(x,A,b):
    return A*np.exp(-b*x)+1

def mytheory(d,delta):
    return 1-np.pi*delta*d**(-2)

# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
rsphere  = 1.0 # 0.8 # 
#sigma    = 5.0 # Change this for the modelling part
Nintervals = 10

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
old_bulk      = False

endlocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/Nocut/'
if moresigmas==True:
    endlocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Moresigmas/Nocut/'

basepath_base      = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_staticbrush/'
endlocation_static = basepath_base+'D_vs_d/Nocut/'

bulklocation = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/'

bulkfilename  = bulklocation + 'diffusion_bulk'+filestext
#if bulk_cut==True:
#    bulkfilename = bulkfilename +'_cut.txt'
#else:
#    bulkfilename = bulkfilename +'_uncut.txt'
bulkfilename = bulkfilename +'_new.txt'

## Files to read
# Folders
endlocation_f   = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_forest/D_vs_d/Nocut/'
brushfilename_f = endlocation_f + 'D_vs_d_forest_better_rms_Nestimates%i.txt' % Nintervals

plotname = endlocation + 'MDstraight_and_RWfixedgeom_norefl_rsphere'+str(rsphere)+'.png'
plotname_wm = endlocation + 'MDstraight_and_RWfixedgeom_norefl_rsphere'+str(rsphere)+'_withmodel.png'

#######

# Forest sims:
brushfile_f = open(brushfilename_f, 'r')

lines = brushfile_f.readlines()
N_f = len(lines)-1

# ds
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
    
    spacings_f[j] = float(words[0])
    DRs_f[j] = float(words[1])
    Dzs_f[j] = float(words[3])
    Dparallel_f[j] = float(words[5])
    # Ds, stdv
    DRs_stdv_f[j] = float(words[2])
    Dzs_stdv_f[j] = float(words[4])
    Dparallel_stdv_f[j] = float(words[6])
    
brushfile_f.close()

## RW part:
Nsteps   = 1000   # Increase #Maybe for other values of hitprob
Nreal    = 10000
folder   = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Randomwalk_fixedgeometry/'

infilename_Dbulk = folder+'2D_bulk_Nsteps%i_Nreal%i_D.txt' %(Nsteps,Nreal)
infile_Ddbulk = open(infilename_Dbulk,'r')
line = infile_Ddbulk.readline()
Dbulk = float(line.split()[0])
Dbulk_rms = float(line.split()[0])
infile_Ddbulk.close()

folder   = folder+'No_reflection/'
infilename_Ddbulk = folder+'D_vs_d_Nsteps%i_Nreal%i_rsphere'%(Nsteps,Nreal)+str(rsphere)+'_norefl.txt' 
infile_Ddbulk = open(infilename_Ddbulk,'r')
lines = infile_Ddbulk.readlines()
print('lines:',lines)

print('infilename_Ddbulk:',infilename_Ddbulk)
dRW  = []
DRW  = []
DRW_rms  = []
DRW_bare = []
for line in lines:
    words = line.split()
    d = float(words[0])
    dRW.append(d)
    Dthis = float(words[1])
    DRW.append(Dthis/Dbulk)
    DRW_bare.append(Dthis)
    rms_bare = float(words[2])
    rms_this = Dthis*np.sqrt((Dbulk_rms/Dbulk)**2+(rms_bare/Dthis)**2)
    DRW_rms.append(rms_this)
infile_Ddbulk.close()

dRW  = np.array(dRW)
DRW  = np.array(DRW)
DRW_rms  = np.array(DRW_rms)

###
#Bulk:

bulkfile  = open(bulkfilename, 'r')
# D_R2  sigmaD_R2 b_R2 sigmab_R2; D_z2  sigmaD_z2  b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2
if old_bulk==True: # Worse rms
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
else:
    bulklines = bulkfile.readlines()
    bulkline  = bulklines[1]
    words     = bulkline.split()
    
    # Ds
    DRs_bulk = float(words[1])
    Dzs_bulk = float(words[3])
    Dparallel_bulk = float(words[5])
    
    # Ds, stdv
    DRs_stdv_bulk = float(words[2])
    Dzs_stdv_bulk = float(words[4])
    Dparallel_stdv_bulk = float(words[6])
    
bulkfile.close()

# Divide by bulk:
#### Forest:
for i in range(N_f):
    DRnew = DRs_f[i]/DRs_bulk
    Dznew = Dzs_f[i]/DRs_bulk
    Dparnew = Dparallel_f[i]/DRs_bulk
    # Ds, stdv
    DRs_stdv_f[i] = abs(DRnew)*np.sqrt((DRs_stdv_f[i]/DRs_f[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dzs_stdv_f[i] = abs(Dznew)*np.sqrt((Dzs_stdv_f[i]/Dzs_f[i])**2+(DRs_stdv_bulk/DRs_bulk)**2)
    Dparallel_stdv_f[i] = abs(Dparnew)*np.sqrt((Dparallel_stdv_f[i]/Dparallel_f[i])**2 +(DRs_stdv_bulk/DRs_bulk)**2)
    # Ds
    DRs_f[i] = DRnew
    Dzs_f[i] = Dznew
    Dparallel_f[i] = Dparnew

#print('DRW_bare:',DRW_bare)
#print('Dbulk:',Dbulk)

plt.figure(figsize=(6,5))
plt.fill_between(spacings_f, Dzs_f+Dzs_stdv_f, Dzs_f-Dzs_stdv_f,color='r',alpha=0.2)
plt.fill_between(spacings_f, Dparallel_f+Dparallel_stdv_f, Dparallel_f-Dparallel_stdv_f,color='k',alpha=0.2)
plt.fill_between(dRW, DRW+DRW_rms, DRW-DRW_rms,color='darkslateblue',alpha=0.2) #CC
plt.plot(spacings_f, Dzs_f, '-v',color='r',label=r'$D_{\perp, \mathregular{straight}}$')
plt.plot(spacings_f, Dparallel_f, '-1',color='k',label=r'$D_{\parallel, \mathregular{straight}}$')
plt.plot(dRW, DRW, '-P',color='darkslateblue',label=r'$D_{\mathregular{RW}}$')
plt.xlabel('$d$ (nm)')
plt.ylabel(r'$D/D_{\mathregular{bulk}}$')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(plotname)

spacings_attempt = spacings_f[1:]
Nattempt = len(spacings_f[1:]) # Or 100, maybe? Need to add another line, in that case.
attempt = np.zeros(Nattempt)
for i in range(Nattempt):
    attempt[i] = 1-np.pi*(1/spacings_f[i+1])**2 # Skipping 1

# Perp
coeffs, covs = curve_fit(mytheory,spacings_attempt, Dzs_f[1:]) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
delta_fz = coeffs[0]
rms_delta_fz = np.sqrt(covs[0,0])

theoryfit_fz = mytheory(spacings_attempt,delta_fz)

# Parallel
coeffs, covs = curve_fit(mytheory,spacings_attempt, Dparallel_f[1:]) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
delta_fparallel = coeffs[0]
rms_delta_fparallel = np.sqrt(covs[0,0])

theoryfit_fparallel = mytheory(spacings_attempt,delta_fparallel)

plt.figure(figsize=(6,5))
plt.fill_between(spacings_f, Dzs_f+Dzs_stdv_f, Dzs_f-Dzs_stdv_f,color='r',alpha=0.2)
plt.fill_between(spacings_f, Dparallel_f+Dparallel_stdv_f, Dparallel_f-Dparallel_stdv_f,color='k',alpha=0.2)
plt.fill_between(dRW, DRW+DRW_rms, DRW-DRW_rms,color='darkslateblue',alpha=0.2) #CC
plt.plot(spacings_f[1:], attempt, '-o' ,label=r'$1-\pi \left(\frac{\sigma}{d}\right)^{2}$')
plt.plot(spacings_attempt, theoryfit_fz, '-o' ,label=r'$1-\pi \delta_\perp\left(\frac{\sigma}{d}\right)^{2}$')
plt.plot(spacings_attempt, theoryfit_fparallel, '-o' ,label=r'$1-\pi \delta_\parallel\left(\frac{\sigma}{d}\right)^{2}$')
plt.plot(spacings_f, Dzs_f, '-v',color='r',label=r'$D_{\perp, \mathregular{straight}}$')
plt.plot(spacings_f, Dparallel_f, '-1',color='k',label=r'$D_{\parallel, \mathregular{straight}}$')
plt.plot(dRW, DRW, '-P',color='darkslateblue',label=r'$D_{\mathregular{RW}}$')
plt.xlabel('$d$ (nm)')
plt.ylabel(r'$D/D_{\mathregular{bulk}}$')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(plotname_wm)

#### Plotting Dperp/Dparallel:

spacings_fraction = []
Dzs_fraction      = []
Dzs_diff          = []
Dzs_rms_fraction  = []
Dzs_rms_diff      = []
ones              = []

spacings_few = []
Dzs_fraction_few = []

i = 0
spacing = spacings_f[i]
while spacing<=25:
    spacings_fraction.append(spacing)    
    Dznew = Dzs_f[i]/Dparallel_f[i]
    # Ds, stdv
    Dzs_stdv = abs(Dznew)*np.sqrt((Dzs_stdv_f[i]/Dzs_f[i])**2+(Dparallel_stdv_f[i]/Dparallel_f[i])**2)
    Dzs_stdv_diff = np.sqrt((Dzs_stdv_f[i])**2+(Dparallel_stdv_f[i])**2)
    Dzs_fraction.append(Dznew)
    Dzs_rms_fraction.append(Dzs_stdv)
    Dzs_diff.append(Dzs_f[i]-Dparallel_f[i])
    Dzs_rms_diff.append(Dzs_stdv_diff)
    ones.append(1)
    if spacing>3:
        spacings_few.append(spacing)
        Dzs_fraction_few.append(Dznew)
    i+=1
    spacing = spacings_f[i]

ones = np.array(ones)
Dzs_diff = np.array(Dzs_diff)
Dzs_fraction = np.array(Dzs_fraction)
Dzs_rms_diff = np.array(Dzs_rms_diff)
spacings_few = np.array(spacings_few)
Dzs_fraction_few = np.array(Dzs_fraction_few)
Dzs_rms_fraction = np.array(Dzs_rms_fraction)
spacings_fraction = np.array(spacings_fraction)

coeffs, covs = curve_fit(myexponential,spacings_few, Dzs_fraction_few) # Using all the data in the fit # Should probably cut the first part, but I don't know how much to include
A_f = coeffs[0]
b_f = coeffs[1]
rms_A_f = np.sqrt(covs[0,0])
rms_b_f = np.sqrt(covs[1,1])

fit_f = myexponential(spacings_few,A_f,b_f)

plt.figure(figsize=(6,5))
plt.fill_between(spacings_fraction, Dzs_fraction+Dzs_rms_fraction, Dzs_fraction-Dzs_rms_fraction,color='lightcoral',alpha=0.2)
plt.plot(spacings_fraction, Dzs_fraction, '-o',color='lightcoral',label='data')
plt.plot(spacings_fraction, ones, '--',color='lightcoral')
plt.plot(spacings_few,fit_f,'--',label=r'%.1f $e^{-%.2fd}$+1' % (A_f,b_f))
plt.xlabel(r'$d$ (nm)')
plt.ylabel(r'$D_{\perp,straight}/D_{\parallel,straight}$')
plt.legend(loc='upper right')

plt.figure(figsize=(6,5))
plt.fill_between(spacings_fraction, Dzs_diff+Dzs_rms_diff, Dzs_diff-Dzs_rms_diff,color='lightcoral',alpha=0.2)
plt.plot(spacings_fraction, Dzs_diff, '-o',color='lightcoral')
plt.xlabel(r'$d$ (nm)')
plt.ylabel(r'$D_{\perp,,straight}/D_{bulk}-D_{\parallel,straight}/D_{bulk}$')
####


print('delta_fz:',delta_fz, '; rms_delta_fz:',rms_delta_fz)
print('delta_fparallel:',delta_fparallel, '; rms_delta_fparallel:',rms_delta_fparallel)

plt.show()