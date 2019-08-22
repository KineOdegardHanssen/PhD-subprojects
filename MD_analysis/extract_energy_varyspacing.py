import matplotlib.pyplot as plt               # To plot
from scipy.optimize import curve_fit
import numpy as np
import random
import math
import time


def find_rms(inarray):
    rms     = 0
    N       = len(inarray)
    average = np.mean(inarray)
    for k in range(N):
        rms += (inarray[k]-average)**2 
    rms = np.sqrt(rms/N)
    return average, rms

start_time_first = time.process_time()

spacings = [1, 2, 3, 4, 5, 10, 40, 100]

all_steps      = []
all_times      = []
all_potEns     = []
all_kinEns     = []
all_totEns     = []
all_bendEns    = []
all_stretchEns = []

etots_av  = np.zeros(len(spacings))
etots_rms = np.zeros(len(spacings))

#plotname     = 'totenergy_harmonic_bond_and_angle_Langevin_nowall_K500_theta0is180_Kbond2000_dtfacs.png'
#plotname2    = 'totenergy_av_rms_harmonic_bond_and_angle_Langevin_nowall_K500_theta0is180_Kbond2000_dtfacs.png'
efilename    = 'energydata_chaingrid_quadratic_M9N101_gridspacing%i_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T310_theta0is180_twofirst_are_fixed_varyspacings.txt'
#plotbasname  = 'basenergies_harmonic_bond_and_angle_Langevin_nowall_K500_theta0is180_Kbond2000_dtfacs.png'
#plotdiffname = 'energydiff_harmonic_bond_and_angle_Langevin_nowall_K500_theta0is180_Kbond2000_dtfacs.png'
efile = open(efilename, 'w') 

for j in range(len(spacings)):
    spacing     = spacings[j] 
    infilename  = 'log.chaingrid_quadratic_M9N101_gridspacing%i_Langevin_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T310_theta0is180_twofirst_are_fixed' % spacing
    
    # The log file in LAMMPS contains a lot of lines detailing the initiation of the system. We want to access the information that comes after that.
    linestart           = 99 #17                  # The line the data starts at.
    printeverynthstep   = 100
    startat             = 50                      # To equilibrate. The number of ns we want to skip

    
    #### Automatic part
    infile = open(infilename, "r")
    lines = infile.readlines()
    # Getting the number of lines, etc.
    totlines = len(lines)         # Total number of lines
    lineend = totlines-1          # Index of last element
    
    steps     = []
    times     = []
    temps     = []
    potEn     = []
    kinEn     = []
    totEn     = []
    bendEn    = []
    stretchEn = []
    coulEn    = []
    
    etot_av = 0
    
    skipelem = 0#10000#10000#90000
    #totlines = N+9+9
    # Extracting data
    i = linestart
    while i<totlines:
        words = lines[i].split()
        if words[0]!='Loop' and words[0]!='WARNING:':
            # Find properties
            #print(words)
            print(words[0])
            steps.append(float(words[0]))
            times.append(float(words[1]))
            temps.append(float(words[2]))
            potEn.append(float(words[3]))
            kinEn.append(float(words[4]))
            totEn.append(float(words[5]))
            bendEn.append(float(words[6]))
            stretchEn.append(float(words[7]))
            coulEn.append(float(words[8]))
            etot_av += totEn[-1]
            i+=1
        else:
            i = totlines # Just in case
            break
    
    infile.close()
    Nt = len(totEn)
    etot_av /= Nt
    
    etot_av, etot_rms           = find_rms(totEn)
    potEn_av, potEn_rms         = find_rms(potEn)
    kinEn_av, kinEn_rms         = find_rms(kinEn)
    bendEn_av, bendEn_rms       = find_rms(bendEn)
    stretchEn_av, stretchEn_rms = find_rms(stretchEn)
    coulEn_av, coulEn_rms       = find_rms(coulEn)
    
    etots_av[j]  = etot_av
    etots_rms[j] = etot_rms
    #					    1	  2	   3	    4		5	6	7	8	   9	      10	11		12	13	
    #		    1     2    3     4      5     6     7     8     9    10    11    12    13
    efile.write('%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n' % (spacing,etot_av,etot_rms, potEn_av, potEn_rms, kinEn_av, kinEn_rms, bendEn_av, bendEn_rms, stretchEn_av, stretchEn_rms, coulEn_av, coulEn_rms))
    
    all_steps.append(steps)
    all_times.append(times)
    all_potEns.append(potEn)
    all_kinEns.append(kinEn)
    all_totEns.append(totEn)
    all_bendEns.append(bendEn)
    all_stretchEns.append(stretchEn)

efile.close()    
'''    
plt.figure(figsize=(6,5))
for i in range(len(dtfacs)):
    times = np.array(all_times[i])
    engys = np.array(all_totEns[i])
    plt.plot(times, engys, label='dt=%.1e' % dts[i])
plt.xlabel(r'Time [ns]', fontsize=16)
plt.ylabel(r'Energy [$10^{-21}$J]', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
plt.legend(loc="upper right")
plt.title(r'Total energy vs time', fontsize=16)
plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.errorbar(dts, etots_av, yerr=etots_rms)
plt.xlabel(r'Time step dt [ns]', fontsize=16)
plt.ylabel(r'Average energy [$10^{-21}$J]', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
plt.title(r'Average energy vs time step', fontsize=16)
plt.savefig(plotname2)

plt.figure(figsize=(6,5))
plt.plot(times, bendEn, label='Bending energy')
plt.plot(times, stretchEn, label='Stretching energy')
plt.xlabel(r'Time [ns]', fontsize=16)
plt.ylabel(r'Energy [$10^{-21}$J]', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
plt.legend(loc="upper right")
plt.title(r'Energy vs time', fontsize=16)
plt.savefig(plotbasname)

plt.figure(figsize=(6,5))
plt.plot(times, potEn-(bendEn+stretchEn))#, label='Bending energy')
#plt.plot(times, stretchEn, label='Stretching energy')
plt.xlabel(r'Time [ns]', fontsize=16)
plt.ylabel(r'Energy [$10^{-21}$J]', fontsize=16)
plt.tight_layout(pad=2.0)#, w_pad=0.0, h_pad=0.5)
#plt.legend(loc="upper right")
plt.title(r'$E_{pot}-E_{bend}-E_{stretch}$ vs time', fontsize=16)
plt.savefig(plotdiffname)
'''
