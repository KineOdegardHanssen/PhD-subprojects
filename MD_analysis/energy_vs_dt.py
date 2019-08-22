import matplotlib.pyplot as plt               # To plot
from scipy.optimize import curve_fit
import numpy as np
import random
import math
import time

start_time_first = time.process_time()

dtfacs = [1.0, 10.0, 100.0]#, 1000.0]
dts    = np.zeros(len(dtfacs))

all_steps      = []
all_times      = []
all_potEns     = []
all_kinEns     = []
all_totEns     = []
all_bendEns    = []
all_stretchEns = []

etots_av  = np.zeros(len(dtfacs))
etots_rms = np.zeros(len(dtfacs))

plotname     = 'totenergy_harmonic_bond_and_angle_Langevin_nowall_K500_theta0is180_Kbond2000_dtfacs.png'
plotname2    = 'totenergy_av_rms_harmonic_bond_and_angle_Langevin_nowall_K500_theta0is180_Kbond2000_dtfacs.png'
efilename    = 'energydata_harmonic_bond_and_angle_Langevin_nowall_K500_theta0is180_Kbond2000_dtfacs.txt'
#plotbasname  = 'basenergies_harmonic_bond_and_angle_Langevin_nowall_K500_theta0is180_Kbond2000_dtfacs.png'
#plotdiffname = 'energydiff_harmonic_bond_and_angle_Langevin_nowall_K500_theta0is180_Kbond2000_dtfacs.png'
efile = open(efilename, 'w') 

for j in range(len(dtfacs)):
    infilename   = 'log.harmonic_bond_and_angle_Langevin_nowall_K500_theta0is180_Kbond2000_dtfac%.1f_atomstylefull' % dtfacs[j]
    
    # The log file in LAMMPS contains a lot of lines detailing the initiation of the system. We want to access the information that comes after that.
    linestart           = 24 #17                  # The line the data starts at.
    printeverynthstep   = 100
    startat             = 50                      # To equilibrate. The number of ns we want to skip
    dt                  = 0.00045/dtfacs[j]       # The time step of our simulation. 0.00045 ns default for nano
    dts[j]              = dt
    
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
    
    etot_av = 0
    
    skipelem = 0#10000#10000#90000
    #totlines = N+9+9
    # Extracting data
    i = linestart
    while i<totlines:
        words = lines[i].split()
        if words[0]!='Loop':
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
            etot_av += totEn[-1]
            i+=1
        else:
            i = totlines # Just in case
            break
    
    infile.close()
    Nt = len(totEn)
    etot_av /= Nt
    
    etot_rms = 0
    for k in range(Nt):
        etot_rms += (totEn[k]-etot_av)**2 
    etot_rms = np.sqrt(etot_rms/Nt)
    
    etots_av[j]  = etot_av
    etots_rms[j] = etot_rms
    
    efile.write('%.16e %.16e %.16e\n' % (dt,etot_av,etot_rms))
    
    all_steps.append(steps)
    all_times.append(times)
    all_potEns.append(potEn)
    all_kinEns.append(kinEn)
    all_totEns.append(totEn)
    all_bendEns.append(bendEn)
    all_stretchEns.append(stretchEn)
    
    
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

'''
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
