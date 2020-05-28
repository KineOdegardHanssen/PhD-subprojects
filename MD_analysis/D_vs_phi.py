from pylab import *
import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
import numpy as np
import math


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

# Functions for curve_fit
def empirical_sand(phi,Dbulk, A,m,n):
    return Dbulk*(A*phi**(1-m))**n

def empirical_soils(phi,Dbulk, B):
    return Dbulk/(phi+B*(1-phi))

def empirical_sediments(phi,Dbulk, C):
    return Dbulk/(1-C*log(phi))

def Deff_theory(Dbulk, phi, modelname):
    ## Models as listed in Shen and Chen, 2007
    #  Models involving phi directly:
    if modelname=='Sherwood':
        print('In %s' % modelname)
        return phi**(3./2)*Dbulk
    if modelname=='PopovicaBrusseau':
        print('In %s' % modelname)
        return phi**(2./3)*Dbulk
    # D/tau**2-models:
    if modelname=='packings_tau2':
        print('In %s' % modelname)
        return Dbulk/(2-phi)
    if modelname=='notmonosized_tau2':
        print('In %s' % modelname)
        return Dbulk*phi**(1./2)
    if modelname=='partsaturated_tau2':
        print('In %s' % modelname)
        return Dbulk*phi**(1./3)
    if modelname=='overlappingspheres_tau2':
        print('In %s' % modelname)
        return Dbulk/(1-log(phi/2))
    if modelname=='overlapping_cylinders_tau2':
        print('In %s' % modelname)
        return Dbulk/(1-log(phi))
    if modelname=='catalyst_tau2':
        print('In %s' % modelname)
        return Dbulk*(1-(1-phi)**(1./3))/phi   
    if modelname=='cationexchange_tau2':
        print('In %s' % modelname)
        return Dbulk*phi**2/(2-phi)**2

psigma = 1

## Filenames
# In-files
basepath     = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/'
mainlocation = basepath+'Diffusion_bead_near_grid/D_vs_d/Brush/Sigma_bead_' +str(psigma) + '/Nocut/'
Dfilename    = mainlocation+'D_vs_d.txt'
phifilen_ws  = basepath + 'd_vs_phi_wsubstrate.txt'
phifilen_oc  = basepath + 'd_vs_phi_withoutsubstrate.txt' 

# Out-files
outfile               = mainlocation + 'D_vs_phi.txt'
plotname_block_oc     = mainlocation + 'D_vs_phi_dyn_nocut_donlychains_block'
plotname_bead_oc      = mainlocation + 'D_vs_phi_dyn_nocut_donlychains_bead'
plotname_block_ws     = mainlocation + 'D_vs_phi_dyn_nocut_withsubstrate_block'
plotname_bead_ws      = mainlocation + 'D_vs_phi_dyn_nocut_withsubstrate_bead'
plotname_block_oc_emp = mainlocation + 'D_vs_phi_dyn_nocut_donlychains_block_empirical'
plotname_bead_oc_emp  = mainlocation + 'D_vs_phi_dyn_nocut_donlychains_bead_empirical'
plotname_block_ws_emp = mainlocation + 'D_vs_phi_dyn_nocut_withsubstrate_block_empirical'
plotname_bead_ws_emp  = mainlocation + 'D_vs_phi_dyn_nocut_withsubstrate_bead_empirical'
plotname_models       = mainlocation + 'D_vs_phi_models'
parameters_empirical  = mainlocation + 'parameters_empirical.txt'

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


## Picking the model
Dbulk = DRs[-1]
# Names
modelname='Sherwood'                     # Bad fit
#modelname='PopovicaBrusseau'             # Bad fit   
#modelname='packings_tau2'                # Bad fit   
#modelname='notmonosized_tau2'            # Bad fit
#modelname='partsaturated_tau2'           # Bad fit
#modelname='overlappingspheres_tau2'      # Bad fit
#modelname='overlapping_cylinders_tau2'   # Bad fit
#modelname='catalyst_tau2'                # At least this one visibly curves
#modelname='cationexchange_tau2'          # Bad fit
theory_block_oc = Deff_theory(Dbulk, phi_block_oc, modelname)
theory_bead_oc = Deff_theory(Dbulk, phi_bead_oc, modelname)
theory_block_ws = Deff_theory(Dbulk, phi_block_ws, modelname)
theory_bead_ws = Deff_theory(Dbulk, phi_bead_ws, modelname)

# Empirical model

#empirical_sand(phi,Dbulk, A,m,n)
#empirical_soils(phi,Dbulk, B)
#empirical_sediments(phi,Dbulk, C)
# Sand
popt, pcov           = curve_fit(lambda phi, A, m, n: empirical_sand(phi,Dbulk, A,m,n), phi_block_oc, DRs)#, maxfev=thismaxfev)                                                                
A_sand_block_oc      = popt[0]
m_sand_block_oc      = popt[1]
n_sand_block_oc      = popt[2]
stdv_A_sand_block_oc = np.sqrt(pcov[0])
stdv_m_sand_block_oc = np.sqrt(pcov[1])
stdv_n_sand_block_oc = np.sqrt(pcov[2])
##
popt, pcov          = curve_fit(lambda phi, A, m, n: empirical_sand(phi,Dbulk, A,m,n), phi_bead_oc, DRs)#, maxfev=thismaxfev)                                                                
A_sand_bead_oc      = popt[0]
m_sand_bead_oc      = popt[1]
n_sand_bead_oc      = popt[2]
stdv_A_sand_bead_oc = np.sqrt(pcov[0])
stdv_m_sand_bead_oc = np.sqrt(pcov[1])
stdv_n_sand_bead_oc = np.sqrt(pcov[2])
##
popt, pcov           = curve_fit(lambda phi, A, m, n: empirical_sand(phi,Dbulk, A,m,n), phi_block_ws, DRs)#, maxfev=thismaxfev)                                                                
A_sand_block_ws      = popt[0]
m_sand_block_ws      = popt[1]
n_sand_block_ws      = popt[2]
stdv_A_sand_block_ws = np.sqrt(pcov[0])
stdv_m_sand_block_ws = np.sqrt(pcov[1])
stdv_n_sand_block_ws = np.sqrt(pcov[2])
##
popt, pcov          = curve_fit(lambda phi, A, m, n: empirical_sand(phi,Dbulk, A,m,n), phi_bead_ws, DRs)#, maxfev=thismaxfev)                                                                
A_sand_bead_ws      = popt[0]
m_sand_bead_ws      = popt[1]
n_sand_bead_ws      = popt[2]
stdv_A_sand_bead_ws = np.sqrt(pcov[0])
stdv_m_sand_bead_ws = np.sqrt(pcov[1])
stdv_n_sand_bead_ws = np.sqrt(pcov[2])


#Soils
popt, pcov           = curve_fit(lambda phi, B: empirical_soils(phi,Dbulk,B), phi_block_oc, DRs)#, maxfev=thismaxfev)                                                                
B_soil_block_oc      = popt[0]
stdv_B_soil_block_oc = np.sqrt(pcov[0])
##
popt, pcov           = curve_fit(lambda phi, B: empirical_soils(phi,Dbulk,B), phi_bead_oc, DRs)#, maxfev=thismaxfev)                                                                
B_soil_bead_oc       = popt[0]
stdv_B_soil_bead_oc  = np.sqrt(pcov[0])
##
popt, pcov           = curve_fit(lambda phi, B: empirical_soils(phi,Dbulk,B), phi_block_ws, DRs)#, maxfev=thismaxfev)                                                                
B_soil_block_ws      = popt[0]
stdv_B_soil_block_ws = np.sqrt(pcov[0])
##
popt, pcov           = curve_fit(lambda phi, B: empirical_soils(phi,Dbulk,B), phi_bead_ws, DRs)#, maxfev=thismaxfev)                                                                
B_soil_bead_ws       = popt[0]
stdv_B_soil_bead_ws  = np.sqrt(pcov[0])


# Sediments
popt, pcov                = curve_fit(lambda phi, C: empirical_sediments(phi,Dbulk,C), phi_block_oc, DRs)#, maxfev=thismaxfev)                                                                
C_sediments_block_oc      = popt[0]
stdv_C_sediments_block_oc = np.sqrt(pcov[0])
##
popt, pcov                = curve_fit(lambda phi, C: empirical_sediments(phi,Dbulk,C), phi_bead_oc, DRs)#, maxfev=thismaxfev)                                                                
C_sediments_bead_oc       = popt[0]
stdv_C_sediments_bead_oc  = np.sqrt(pcov[0])
##
popt, pcov                = curve_fit(lambda phi, C: empirical_sediments(phi,Dbulk,C), phi_block_ws, DRs)#, maxfev=thismaxfev)                                                                
C_sediments_block_ws      = popt[0]
stdv_C_sediments_block_ws = np.sqrt(pcov[0])
##
popt, pcov                = curve_fit(lambda phi, C: empirical_sediments(phi,Dbulk,C), phi_bead_ws, DRs)#, maxfev=thismaxfev)                                                                
C_sediments_bead_ws       = popt[0]
stdv_C_sediments_bead_ws  = np.sqrt(pcov[0])

# Empirical curves
#Sand
emp_sand_block_oc = empirical_sand(phi_block_oc,Dbulk, A_sand_block_oc,m_sand_block_oc,n_sand_block_oc)
emp_sand_bead_oc = empirical_sand(phi_bead_oc,Dbulk, A_sand_bead_oc,m_sand_bead_oc,n_sand_bead_oc)
emp_sand_block_ws = empirical_sand(phi_block_ws,Dbulk,A_sand_block_ws,m_sand_block_ws,n_sand_block_ws)
emp_sand_bead_ws = empirical_sand(phi_bead_ws,Dbulk, A_sand_bead_ws,m_sand_bead_ws,n_sand_bead_ws)
#Soil
emp_soil_block_oc    = empirical_soils(phi_block_oc,Dbulk, B_soil_block_oc)
emp_soil_bead_oc     = empirical_soils(phi_bead_oc,Dbulk, B_soil_bead_oc)
emp_soil_block_ws    = empirical_soils(phi_block_ws,Dbulk,B_soil_block_ws)
emp_soil_bead_ws     = empirical_soils(phi_bead_ws,Dbulk, B_soil_bead_ws)
#Sediments
emp_sediments_block_oc    = empirical_sediments(phi_block_oc,Dbulk, C_sediments_block_oc)
emp_sediments_bead_oc     = empirical_sediments(phi_bead_oc,Dbulk, C_sediments_bead_oc)
emp_sediments_block_ws    = empirical_sediments(phi_block_ws,Dbulk,C_sediments_block_ws)
emp_sediments_bead_ws     = empirical_sediments(phi_bead_ws,Dbulk, C_sediments_bead_ws)
# Deviations from data:
# Sand
emp_sand_block_oc_rmsd = rmsd(emp_sand_block_oc,DRs)
emp_sand_bead_oc_rmsd  = rmsd(emp_sand_bead_oc,DRs)
emp_sand_block_ws_rmsd = rmsd(emp_sand_block_ws,DRs)
emp_sand_bead_ws_rmsd  = rmsd(emp_sand_bead_ws,DRs)
# Soil
emp_soil_block_oc_rmsd = rmsd(emp_soil_block_oc,DRs)
emp_soil_bead_oc_rmsd  = rmsd(emp_soil_bead_oc,DRs)
emp_soil_block_ws_rmsd = rmsd(emp_soil_block_ws,DRs)
emp_soil_bead_ws_rmsd  = rmsd(emp_soil_bead_ws,DRs)
# Sediments
emp_sediments_block_oc_rmsd = rmsd(emp_sediments_block_oc,DRs)
emp_sediments_bead_oc_rmsd  = rmsd(emp_sediments_bead_oc,DRs)
emp_sediments_block_ws_rmsd = rmsd(emp_sediments_block_ws,DRs)
emp_sediments_bead_ws_rmsd  = rmsd(emp_sediments_bead_ws,DRs)

# Write params and rmsd to file
outfile_paramemp = open(parameters_empirical,'w')
outfile_paramemp.write('Sand (A,m,n) --> D/tau^2, tau^2 = (A*phi^(1-m)^n):\n')
outfile_paramemp.write('Block, only chains: %.16e %.16e %.16e rmsd %.16e\n' % (A_sand_block_oc,m_sand_block_oc,n_sand_block_oc,emp_sand_block_oc_rmsd))
outfile_paramemp.write('Bead, only chains: %.16e %.16e %.16e rmsd %.16e\n' % (A_sand_bead_oc,m_sand_bead_oc,n_sand_bead_oc,emp_sand_bead_oc_rmsd))
outfile_paramemp.write('Block, with substrate: %.16e %.16e %.16e rmsd %.16e\n' % (A_sand_block_ws,m_sand_block_ws,n_sand_block_ws,emp_sand_block_ws_rmsd))
outfile_paramemp.write('Bead, with substrate: %.16e %.16e %.16e rmsd %.16e\n' % (A_sand_bead_ws,m_sand_bead_ws,n_sand_bead_ws,emp_sand_bead_ws_rmsd))
outfile_paramemp.write('\nSoils (B) --> D/tau^2, tau^2 = (phi+B(1-phi)):\n')
outfile_paramemp.write('Block, only chains: %.16e rmsd %.16e\n' % (B_soil_block_oc,emp_soil_block_oc_rmsd))
outfile_paramemp.write('Bead, only chains: %.16e rmsd %.16e\n' % (B_soil_bead_oc,emp_soil_bead_oc_rmsd))
outfile_paramemp.write('Block, with substrate: %.16e rmsd %.16e\n' % (B_soil_block_ws,emp_soil_block_ws_rmsd))
outfile_paramemp.write('Bead, with substrate: %.16e rmsd %.16e\n' % (B_soil_bead_ws,emp_soil_bead_ws_rmsd))
outfile_paramemp.write('\nSediments (C) --> D/tau^2, tau^2 = 1-C*ln(phi):\n')
outfile_paramemp.write('Block, only chains: %.16e rmsd %.16e\n' % (C_sediments_block_oc,emp_sediments_block_oc_rmsd))
outfile_paramemp.write('Bead, only chains: %.16e rmsd %.16e\n' % (C_sediments_bead_oc,emp_sediments_bead_oc_rmsd))
outfile_paramemp.write('Block, with substrate: %.16e rmsd %.16e\n' % (C_sediments_block_ws,emp_sediments_block_ws_rmsd))
outfile_paramemp.write('Bead, with substrate: %.16e rmsd %.16e\n' % (C_sediments_bead_ws,emp_sediments_bead_ws_rmsd))


# Plot of empirical vs data
# Only chains in phi
plt.figure(figsize=(10,5))
ax = plt.subplot(111)
ax.plot(phi_block_oc, DRs, label=r'Data')
ax.plot(phi_block_oc, emp_sand_block_oc, '--', label=r'Sand, empirical')
ax.plot(phi_block_oc, emp_soil_block_oc, '--', label=r'Soil, empirical')
ax.plot(phi_block_oc, emp_sediments_block_oc, '--', label=r'Sediment, empirical')
plt.xlabel(r'Porosity $\phi$ (block, no substr.)')
plt.ylabel(r'Diffusion constant $D$')
plt.title(r'$D$ vs $\phi$, block, no substr.')
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.axis([0, xmax_plot, 0, max(dR2[0:xmax_plot])])
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_block_oc_emp)

plt.figure(figsize=(10,5))
ax = plt.subplot(111)
ax.plot(phi_bead_oc, DRs, label=r'Data')
ax.plot(phi_bead_oc, emp_sand_bead_oc, '--', label=r'Sand, empirical')
ax.plot(phi_bead_oc, emp_soil_bead_oc, '--', label=r'Soil, empirical')
ax.plot(phi_bead_oc, emp_sediments_bead_oc, '--', label=r'Sediment, empirical')
plt.xlabel(r'Porosity $\phi$ (bead, no substr.)')
plt.ylabel(r'Diffusion constant $D$')
plt.title(r'$D$ vs $\phi$, bead, no substr.')
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_bead_oc_emp)

# With substrate
plt.figure(figsize=(10,5))
ax = plt.subplot(111)
ax.plot(phi_block_ws, DRs, label=r'Data')
ax.plot(phi_block_ws, emp_sand_block_ws, '--', label=r'Sand, empirical')
ax.plot(phi_block_ws, emp_soil_block_ws, '--', label=r'Soil, empirical')
ax.plot(phi_block_ws, emp_sediments_block_ws, '--', label=r'Sediment, empirical')
plt.xlabel(r'Porosity $\phi$ (block, substr.)')
plt.ylabel(r'Diffusion constant $D$')
plt.title(r'$D$ vs $\phi$, block, substr.')
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_block_ws_emp)

plt.figure(figsize=(10,5))
ax = plt.subplot(111)
ax.plot(phi_bead_ws, DRs, label=r'Data')
ax.plot(phi_bead_ws, emp_sand_bead_ws, '--', label=r'Sand, empirical')
ax.plot(phi_bead_ws, emp_soil_bead_ws, '--', label=r'Soil, empirical')
ax.plot(phi_bead_ws, emp_sediments_bead_ws, '--', label=r'Sediment, empirical')
plt.xlabel(r'Porosity $\phi$ (bead, substr.)')
plt.ylabel(r'Diffusion constant $D$')
plt.title(r'$D$ vs $\phi$, bead, substr.')
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_bead_ws_emp)


# Plotting models
longphis = np.linspace(0.1,1,501)
theory_Sherwood                   = Deff_theory(Dbulk, longphis, 'Sherwood')                    # 1
theory_PopovicaBrusseau           = Deff_theory(Dbulk, longphis, 'PopovicaBrusseau')            # 2
theory_packings_tau2              = Deff_theory(Dbulk, longphis, 'packings_tau2')               # 3
theory_notmonosized_tau2          = Deff_theory(Dbulk, longphis, 'notmonosized_tau2')           # 4
theory_partsaturated_tau2         = Deff_theory(Dbulk, longphis, 'partsaturated_tau2')          # 5
theory_overlappingspheres_tau2    = Deff_theory(Dbulk, longphis, 'overlappingspheres_tau2')     # 6
theory_overlapping_cylinders_tau2 = Deff_theory(Dbulk, longphis, 'overlapping_cylinders_tau2')  # 7
theory_catalyst_tau2              = Deff_theory(Dbulk, longphis, 'catalyst_tau2')               # 8
theory_cationexchange_tau2        = Deff_theory(Dbulk, longphis, 'cationexchange_tau2')         # 9

# Test of theory graphs
plt.figure(figsize=(10,5))
ax = plt.subplot(111)
ax.plot(longphis, theory_Sherwood, label=r'Sherwood')                                # 1
ax.plot(longphis, theory_PopovicaBrusseau, label=r'Popovica Brusseau')               # 2
ax.plot(longphis, theory_packings_tau2, label=r'Isotropic packings')                 # 3
ax.plot(longphis, theory_notmonosized_tau2, label=r'Not monosized spheres')          # 4
ax.plot(longphis, theory_partsaturated_tau2, label=r'Partly saturated')              # 5
ax.plot(longphis, theory_overlappingspheres_tau2, label=r'Overlapping spheres')      # 6
ax.plot(longphis, theory_overlapping_cylinders_tau2, label=r'Overlapping cylinders') # 7
ax.plot(longphis, theory_catalyst_tau2, label=r'Catalyst')                           # 8
ax.plot(longphis, theory_cationexchange_tau2, label=r'Cation exchange')              # 9
plt.xlabel(r'Porosity $\phi$')
plt.ylabel(r'Diffusion constant $D$')
plt.title(r'$D$ vs $\phi$ for different models')
plt.tight_layout()
plt.legend(loc='upper left')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(plotname_models)



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



