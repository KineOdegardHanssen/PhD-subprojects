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

def theory(x,A,b):
    return A*x**b

Nintervals = 10

# Input parameters for file selection: # I will probably add more, but I want to make sure the program is running first
psigmas  = [0.25,0.5,1,1.5,2,3,5,10]#,25]
spacing  = 10
pmass    = 1 # I don't use this anymore, but it got stuck in the file names. Have constant mass density now.
damp     = 10
N        = len(psigmas)

# Input booleans for file selection:
bulkdiffusion = True
substrate     = False
confignrs     = np.arange(1,1001)

# Ds
DRs = np.zeros(N)
Dxs = np.zeros(N)
Dys = np.zeros(N)
Dzs = np.zeros(N)
Dparallel = np.zeros(N)
# Ds, stdv
DRs_stdv = np.zeros(N)
Dxs_stdv = np.zeros(N)
Dys_stdv = np.zeros(N)
Dzs_stdv = np.zeros(N)
Dparallel_stdv = np.zeros(N)

# bs
bRs = np.zeros(N)
bxs = np.zeros(N)
bys = np.zeros(N)
bzs = np.zeros(N)
bparallel = np.zeros(N)
# bs, stdv
bRs_stdv = np.zeros(N)
bxs_stdv = np.zeros(N)
bys_stdv = np.zeros(N)
bzs_stdv = np.zeros(N)
bparallel_stdv = np.zeros(N)


if bulkdiffusion==True:
    parentfolder = 'Bulk/'
    filestext    = '_seed'+str(confignrs[0])+'to'+str(confignrs[-1])
    systemtype   = 'bulk'
    if substrate==True:
        parentfolder = 'Bulk_substrate/'
        systemtype   = 'substrate'
else:
    parentfolder = 'Brush/'
    systemtype   = 'brush'
    filestext    = 'config'+str(confignrs[0])+'to'+str(confignrs[-1])

endlocation_out = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/D_vs_d/'+parentfolder+'Varysigmas/'
outfilename  = endlocation_out+'D_vs_sigma_better_rms_Nestimates%i_pmass1.txt' % Nintervals
plotname     = endlocation_out+'D_vs_sigma_bulk_better_rms_Nestimates%i_pmass1.png' % Nintervals
plotname_two = endlocation_out+'D_vs_sigma_bulk_better_rms_Nestimates%i_twoinone_pmass1.png' % Nintervals
plotname_fit = endlocation_out+'D_vs_sigma_fit_better_rms_Nestimates%i_pmass1.png' % Nintervals
indfilename  = endlocation_out+'D_vs_sigma_fitindices_better_rms_Nestimates%i_pmass1.txt' % Nintervals

outfile = open(outfilename, 'w')
outfile.write('d   D_R2   sigmaD_R2  b_R2 sigmab_R2; D_z2  sigmaD_z2 b_z2  sigmaD_z2; D_par2 sigmaD_par2  b_par2  sigmab_par2\n')

indexfile = open(indfilename, 'w')
indexfile.write('Start_index_R     end_index_R     Start_index_ort     end_index_ort     Start_index_par     end_index_par\n')

for i in range(N):
    psigma = psigmas[i]
    endlocation_in   = '/home/kine/Documents/Backup2_P2_PolymerMD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/' % damp +'Pure_bulk/'+ 'Sigma_bead_' +str(psigma) + '/'
    infilename = endlocation_in+'diffusion'+filestext+'_better_rms_Nestimates%i_pmass1.txt' % Nintervals
    metaname   = endlocation_in+'indices_for_fit.txt'

    #print('infilename_all:',infilename_all)
    
    # Read in:
    #### Automatic part
    ## Find the extent of the polymers: Max z-coord of beads in the chains
    infile = open(infilename, "r")
    lines = infile.readlines() # This takes some time
    # Getting the number of lines, etc.
    line = lines[1]
    words = line.split()
    
    # Ds
    DRs[i] = float(words[0])
    Dzs[i] = float(words[2])
    Dparallel[i] = float(words[4])
    # Ds, stdv
    DRs_stdv[i] = float(words[1])
    Dzs_stdv[i] = float(words[3])
    Dparallel_stdv[i] = float(words[5])
    
    infile.close()
    
    outfile.write('%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n' % (psigma, DRs[i], DRs_stdv[i],  Dzs[i], Dzs_stdv[i], Dparallel[i], Dparallel_stdv[i]))
    
    metafile   = open(metaname, 'r')
    mlines      = metafile.readlines()
    startindex_R  = int(mlines[0].split()[1]) # New setup
    endindex_R    = int(mlines[1].split()[1])
    startindex_z  = int(mlines[2].split()[1])
    endindex_z    = int(mlines[3].split()[1])
    startindex_p  = int(mlines[4].split()[1])
    endindex_p    = int(mlines[5].split()[1])
    metafile.close()
    indexfile.write('%.2f %i %i %i %i %i %i\n' % (spacing, startindex_R, endindex_R, startindex_z, endindex_z, startindex_p, endindex_p))

outfile.close()

psigmas = np.array(psigmas)

# Make fit

# Fit, start
popt, pcov = curve_fit(theory, psigmas[:5], DRs[:5])
A = popt[0]
b = popt[1]
rms_params = np.sqrt(np.diag(pcov))
A_stdv = rms_params[0]
b_stdv = rms_params[1]

theory_graph = theory(psigmas, A, b)

# Fit, end
popt, pcov = curve_fit(theory, psigmas[-3:], DRs[-3:])
Aend = popt[0]
bend = popt[1]
rms_params = np.sqrt(np.diag(pcov))
Aend_stdv = rms_params[0]
bend_stdv = rms_params[1]

theoryend_graph = theory(psigmas, Aend, bend)

# Fit, all:
popt, pcov = curve_fit(theory, psigmas, DRs)
Aall = popt[0]
ball = popt[1]
rms_params = np.sqrt(np.diag(pcov))
Aall_stdv = rms_params[0]
ball_stdv = rms_params[1]

theoryall_graph = theory(psigmas, Aall, ball)
    
plt.figure(figsize=(6,5))
plt.errorbar(psigmas, DRs, yerr=DRs_stdv, capsize=2, label=r'$D_R$')
plt.errorbar(psigmas, Dzs, yerr=Dzs_stdv, capsize=2, label=r'$D_\perp$')
plt.errorbar(psigmas, Dparallel, yerr=Dparallel_stdv, capsize=2, label=r'$D_\parallel$')
plt.xlabel(r'$\sigma_{\mathregular{LJ}}$ (nm)')
plt.ylabel(r'$D$ (m$^2$/s)')
#plt.title('Diffusion constant $D$ vs $\sigma_b$, %s' % systemtype)
plt.tight_layout()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='upper right')
plt.savefig(plotname)

'''
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
ax1.fill_between(psigmas, DRs+DRs_stdv, DRs-DRs_stdv, color='r', alpha=0.3)
ax1.fill_between(psigmas, Dzs+Dzs_stdv, Dzs-Dzs_stdv, color='b', alpha=0.3)
ax1.fill_between(psigmas, Dparallel+Dparallel_stdv, Dparallel-Dparallel_stdv, color='g', alpha=0.3)
ax1.plot(psigmas, DRs, '-s', color='r', label=r'$D_R$')
ax1.plot(psigmas, Dzs, '-o', color='b', label=r'$D_\perp$')
ax1.plot(psigmas, Dparallel, '-*', color='g', label=r'$D_\parallel$')
#ax1.set(xlabel=r'$\sigma_{\mathregular{LJ}} (nm)$',ylabel=r'$D$ (m$^2$/s)')
ax1.set_xlabel(xlabel=r'$\sigma_{\mathregular{LJ}}$ (nm)',labelpad=2)
ax1.set_ylabel(ylabel=r'$D$ (m$^2$/s)')
#fig.suptitle('Diffusion constant $D$ vs $\sigma_b$, %s' % systemtype)
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0),useMathText=True)
ax1.legend(loc='upper right')
ax1.set_title('a)', y=1.05, loc='left')#, pad=8, loc='left')
ax2.loglog(psigmas, DRs, '-s', color='r')#, label=r'$D_R$')
ax2.loglog(psigmas, Dzs, '-o', color='b')#, label=r'$D_\perp$')
ax2.loglog(psigmas, Dparallel, '-*', color='g')#, label=r'$D_\parallel$')
#ax2.loglog(psigmas, theory_graph, '--', label=r'A$\sigma_{\mathregular{LJ}}^{b}$, first points')
#ax2.loglog(psigmas, theoryend_graph, '--', label=r'A$\sigma_{\mathregular{LJ}}^{b}$, end points')
ax2.loglog(psigmas, theoryall_graph, '--', label=r'A$\sigma_{\mathregular{LJ}}^{b}$')
#ax2.set(xlabel=r'$\sigma_{\mathregular{LJ}}$ (nm)',ylabel=r'$D$ (m$^2$/s)')
ax2.set_xlabel(xlabel=r'$\sigma_{\mathregular{LJ}}$ (nm)',labelpad=-1.5)
ax2.set_ylabel(ylabel=r'$D$ (m$^2$/s)',labelpad=-5)
ax2.legend(loc='lower left')
ax2.set_title('b)', y=1.05, loc='left')#, pad=8, loc='left') # [Doesn't work: , set_x=-1]
'''

fig, ax1 = plt.subplots()
ax1.fill_between(psigmas, DRs+DRs_stdv, DRs-DRs_stdv, color='r', alpha=0.3)
ax1.fill_between(psigmas, Dzs+Dzs_stdv, Dzs-Dzs_stdv, color='b', alpha=0.3)
ax1.fill_between(psigmas, Dparallel+Dparallel_stdv, Dparallel-Dparallel_stdv, color='g', alpha=0.3)
ax1.plot(psigmas, DRs, '-s', color='r', label=r'$D_R$')
ax1.plot(psigmas, Dzs, '-o', color='b', label=r'$D_\perp$')
ax1.plot(psigmas, Dparallel, '-*', color='g', label=r'$D_\parallel$')
#ax1.set(xlabel=r'$\sigma_{\mathregular{LJ}} (nm)$',ylabel=r'$D$ (m$^2$/s)')
ax1.set_xlabel(xlabel=r'$\sigma_{\mathregular{LJ}}$ (nm)',labelpad=2)
ax1.set_ylabel(ylabel=r'$D$ (m$^2$/s)')
#fig.suptitle('Diffusion constant $D$ vs $\sigma_b$, %s' % systemtype)
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0),useMathText=True)
ax1.legend(loc='upper right')
#ax1.set_title('a)', y=1.05, loc='left')#, pad=8, loc='left')

plt.tight_layout()

ax2 = fig.add_axes([0.35, 0.35, 0.4, 0.4])
ax2.loglog(psigmas, DRs, '-s', color='r')#, label=r'$D_R$')
ax2.loglog(psigmas, Dzs, '-o', color='b')#, label=r'$D_\perp$')
ax2.loglog(psigmas, Dparallel, '-*', color='g')#, label=r'$D_\parallel$')
#ax2.loglog(psigmas, theory_graph, '--', label=r'A$\sigma_{\mathregular{LJ}}^{b}$, first points')
#ax2.loglog(psigmas, theoryend_graph, '--', label=r'A$\sigma_{\mathregular{LJ}}^{b}$, end points')
ax2.loglog(psigmas, theoryall_graph, '--', label=r'A$\sigma_{\mathregular{LJ}}^{b}$')
#ax2.set(xlabel=r'$\sigma_{\mathregular{LJ}}$ (nm)',ylabel=r'$D$ (m$^2$/s)')
ax2.set_xlabel(xlabel=r'$\sigma_{\mathregular{LJ}}$ (nm)',labelpad=-1.5)
ax2.set_ylabel(ylabel=r'$D$ (m$^2$/s)',labelpad=-5)
ax2.legend(loc='lower left')
#ax2.set_title('b)', y=1.05, loc='left')#, pad=8, loc='left') # [Doesn't work: , set_x=-1]

plt.savefig(plotname_two)

plt.show() # For testing purposes

# Get exponent of smallest d's:
logds  = np.log(psigmas)
logDRs = np.log(DRs)

a = (logDRs[2]-logDRs[0])/(logds[2]-logds[0])
print('Exponent:',a)
print('-------------')
print('Fit performed on: sigma=',psigmas[:5])
print('Exponent, from fit:',b, '+/-', b_stdv)
print('-------------')
print('Fit performed on: sigma=',psigmas[-3:])
print('Exponent, from fit:',bend, '+/-', bend_stdv)
print('-------------')
print('Fit performed on: sigma=',psigmas)
print('Exponent, from fit:',ball, '+/-', ball_stdv)


