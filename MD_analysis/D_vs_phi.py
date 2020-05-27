from pylab import *
import matplotlib.pyplot as plt                     # To plot
import numpy as np
import math


def Deff_theory(Dbulk, phi):
    return phi**(3./2)*Dbulk

psigma = 1

## Filenames
# In-files
basepath     = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/'
mainlocation = basepath+'Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/Nocut/'
Dfilename    = mainlocation+'D_vs_d.txt'
phifilen_ws  = basepath + 'd_vs_phi_wsubstrate.txt'
phifilen_oc  = basepath + 'd_vs_phi_withoutsubstrate.txt' 

# Out-files
outfile = mainlocation + 'D_vs_phi.txt'
plotname_block_oc = mainlocation + 'D_vs_phi_dyn_nocut_donlychains_block'
plotname_bead_oc = mainlocation + 'D_vs_phi_dyn_nocut_donlychains_bead'
plotname_block_ws = mainlocation + 'D_vs_phi_dyn_nocut_withsubstrate_block'
plotname_bead_ws = mainlocation + 'D_vs_phi_dyn_nocut_withsubstrate_bead'

## Read content of files
# Diffusion constant
Dfile = open(Dfilename,'r')
header = Dfile.readline()
lines  = Dfile.readlines()
N      = len(lines)

DRs = np.zeros(N)
#Dxs = np.zeros(N)
#Dys = np.zeros(N)
Dzs = np.zeros(N)
Dparallel = np.zeros(N)
DRs_stdv = np.zeros(N)
Dzs_stdv = np.zeros(N)
Dparallel_stdv = np.zeros(N)

for i in range(N):
    words = lines[i].split()
    
    DRs[i] = float(words[1])
    Dzs[i] = float(words[5])
    Dparallel[i] = float(words[9])
    # Ds, stdv
    DRs_stdv[i] = float(words[2])
    Dzs_stdv[i] = float(words[6])
    Dparallel_stdv[i] = float(words[10])
Dfile.close()

# d (without substrate)
phifile = open(phifilen_oc, 'r')
header  = phifile.readline()
lines   = phifile.readlines()

phi_block_oc = np.zeros(N)
phi_bead_oc  = np.zeros(N)

for i in range(N):
    words = lines[i].split()
    phi_block_oc[i] = float(words[1])
    phi_bead_oc[i]  = float(words[2])
phifile.close()

# d (with substrate)
phifile = open(phifilen_ws, 'r')
header  = phifile.readline()
lines   = phifile.readlines()

phi_block_ws = np.zeros(N)
phi_bead_ws  = np.zeros(N)

for i in range(N):
    words = lines[i].split()
    phi_block_ws[i] = float(words[1])
    phi_bead_ws[i]  = float(words[2])
phifile.close()


#
Dbulk = DRs[-1]
theory_block_oc = Deff_theory(Dbulk, phi_block_oc)
theory_bead_oc = Deff_theory(Dbulk, phi_bead_oc)
theory_block_ws = Deff_theory(Dbulk, phi_block_ws)
theory_bead_ws = Deff_theory(Dbulk, phi_bead_ws)

# Only chains in phi
plt.figure(figsize=(6,5))
plt.plot(phi_block_oc, DRs, label=r'$D$')
plt.plot(phi_block_oc, theory_block_oc, '--', label=r'$D$, theory')
plt.plot(phi_block_oc, Dzs, label=r'$D_\perp$')
plt.plot(phi_block_oc, Dparallel, label=r'$D_\parallel$')
plt.xlabel(r'Porosity $\phi$ (block, no substr.)')
plt.ylabel(r'Diffusion constant $D$')
plt.title(r'$D$ vs $\phi$, block, no substr.')
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.axis([0, xmax_plot, 0, max(dR2[0:xmax_plot])])
plt.savefig(plotname_block_oc)

plt.figure(figsize=(6,5))
plt.plot(phi_bead_oc, DRs, label=r'$D$')
plt.plot(phi_bead_oc, theory_bead_oc, '--', label=r'$D$, theory')
plt.plot(phi_bead_oc, Dzs, label=r'$D_\perp$')
plt.plot(phi_bead_oc, Dparallel, label=r'$D_\parallel$')
plt.xlabel(r'Porosity $\phi$ (bead, no substr.)')
plt.ylabel(r'Diffusion constant $D$')
plt.title(r'$D$ vs $\phi$, bead, no substr.')
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.axis([0, xmax_plot, 0, max(dR2[0:xmax_plot])])
plt.savefig(plotname_bead_oc)

# Substrate and chains in phi
plt.figure(figsize=(6,5))
plt.plot(phi_block_ws, DRs, label=r'$D$')
plt.plot(phi_block_ws, theory_block_ws, '--', label=r'$D$, theory')
plt.plot(phi_block_ws, Dzs, label=r'$D_\perp$')
plt.plot(phi_block_ws, Dparallel, label=r'$D_\parallel$')
plt.xlabel(r'Porosity $\phi$ (block, substr.)')
plt.ylabel(r'Diffusion constant $D$')
plt.title(r'$D$ vs $\phi$, block, substr.')
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.axis([0, xmax_plot, 0, max(dR2[0:xmax_plot])])
plt.savefig(plotname_block_ws)

plt.figure(figsize=(6,5))
plt.plot(phi_bead_ws, DRs, label=r'$D$')
plt.plot(phi_bead_ws, theory_bead_ws, '--', label=r'$D$, theory')
plt.plot(phi_bead_ws, Dzs, label=r'$D_\perp$')
plt.plot(phi_bead_ws, Dparallel, label=r'$D_\parallel$')
plt.xlabel(r'Porosity $\phi$ (bead, substr.)')
plt.ylabel(r'Diffusion constant $D$')
plt.title(r'$D$ vs $\phi$, bead, substr.')
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.axis([0, xmax_plot, 0, max(dR2[0:xmax_plot])])
plt.savefig(plotname_bead_ws)



