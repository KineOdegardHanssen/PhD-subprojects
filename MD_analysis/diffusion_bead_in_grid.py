from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys
import math

# This is a specified script that reads a log file, treats the data and plots the temperature as a function of time. If you want to do something similar with another quantity, please look at basic_reader.py, which is more of a template/sketch of a script.

# This script works when we have not specified the format of the log file from LAMMPS. It could also work when the temperature is given by the fifth column.

#### Input data
filename            = "log.particle_in_chaingrid_quadratic_M9N101_ljunits_gridspacing5_Langevin_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debyecutoff3_chargeelementary-1_effectivedielectric0.00881819074717447_T3_theta0is180_pmass1.5" #"log.diffusion_density0p6_T1p0_right"#"log.diffusion_density0p6_T1p0"#"log.diffusion_density0p8_T0p5"
# The log file in LAMMPS contains a lot of lines detailing the initiation of the system. We want to access the information that comes after that.
linestart           = 244#17            # The line the data starts at.
printeverynthstep   = 0
#startat             = 50   
#skipelem            = 50
dt                  = 0.00045         # The time step of our simulation. 

#### Automatic part
infile = open(filename, "r")
lines = infile.readlines()
# Getting the number of lines, etc.
totlines = len(lines)         # Total number of lines
lineend = totlines-1          # Index of last element


# Readying lists
time    = []
tempr   = [] # Not sure if I need this
msd     = []
nlines = 0

# Extracting data
for i in range(linestart, totlines):
    words = lines[i].split()
    if words[0]=='Loop':
        break
    time.append(float(words[1]))
    tempr.append(float(words[2]))
    msd.append(float(words[10]))
    nlines += 1

print("nlines:", nlines)
# Making the proper arrays
time        = array(time)
tempr       = array(tempr)
msd         = array(msd)
time *= dt

# After equilibration
#time_eq     = time[skipelem:nlines-1]
#tempr_eq    = tempr[skipelem:nlines-1]
#msd_eq      = msd[skipelem:nlines-1]

# Using numpy's polyfit function to fit to line # I hope this works
'''
coeffs = polyfit(time_eq, msd_eq, 1)
a0  = coeffs[0]
D   = a0/6.
# Finding the average temperature (after equilibration)
avtemp = mean(tempr_eq)
'''

coeffs = polyfit(time, msd, 1)
a0  = coeffs[0]
D   = a0/6.
# Finding the average temperature (after equilibration)
avtemp = mean(tempr)


print(" ")
print("Program diffusion_bead_in_grid.py,")
print("dt =", dt)

print(" ")
print("Slope of msd vs t after eq:")
print("slope=", a0)
print("D=", D)

straightline = a0*time+coeffs[1]

Ktemp = avtemp*119.74
print(" ")
print("Input temp=", avtemp, "tau")
print("Temp=",Ktemp, "K")
print("Temp=",Ktemp-273.15, "C")

# Plotting
figure(figsize=(6,5))
plot(time, msd, label='Data')
plot(time, straightline, 'r--', label='Line fit')
xlabel(r'Time [$\tau$]', fontsize=14)
ylabel(r'$<r^2>$ [$\sigma^2$]', fontsize=14)
title(r'$<r^2>$ vs time, T=%.2f $\epsilon/k$, $\rho=0.6$' % avtemp, fontsize=14)
tight_layout()
legend(loc='lower right')
show()
#show(block=False)
