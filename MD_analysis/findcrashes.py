from skimage import morphology, measure, io, util   # For performing morphology operations (should maybe import more from skimage...)
from mpl_toolkits.mplot3d import Axes3D             # Plotting in 3D
import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
from scipy.ndimage import measurements, convolve    # Should have used the command below, but changing now will lead to confusion
from scipy import ndimage                           # For Euclidean distance measurement
import maintools_percolation as perctools
import data_treatment as datr
import numpy as np
import random
import math
import time
import os
import glob
import copy

spacing = 10
sigma_atom = 1
configs = [1000]#np.arange(1,1001)
Nconfigs = len(configs)

basefolder = '/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Brush/Quadr_M9N101_ljunits_Langevin_scaled_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debcutoff3_chargeel-1_effdiel0.00881819074717447_T3_ljcut1p122/Spacing%i/Sigma_bead_'% (spacing)+ str(sigma_atom)+'/'
folder = basefolder + 'log-files/'
outfolder = basefolder
outfilename = outfolder + 'collisiontimes.txt' #_test_every.txt'
outfilename_average = outfolder + 'average_collisioninterval.txt'
all_crashtimes = []

outfile = open(outfilename, 'w')
outfile.write('Results written by findcrashes.py. Order of elements: time of collision, potential energy of bead\n')

skippedfiles = 0
filescounter = 0
av_coll_int  = 0 # Average collision interval
coll_counter = 0 # Total number of collisions
coll_int_every = []
for config in configs: # Write to file inside loop?
    print('config:', config)
    prev_coll = 0 # To find the collision interval
    
    filename_infolder = 'log.confignr%i' % config #+ '_every'
    infilename = folder + filename_infolder
    
    try:
        infile_all = open(infilename, "r")
    except:
        print('Oh, log-file! Where art thou?')
        skippedfiles += 1
        continue # Skipping this file if it does not exist
    # Moving on, if the file
    filescounter += 1
    outfile.write('Config nr: %i\n' % config)
    infile = open(infilename, 'r')
    lines = infile.readlines()
    
    startit = 0
    times   = []
    crashtimes   = []
    pot_energies = []
    for line in lines:
        words = line.split()
        if len(words)>1:
            if words[0]=='Loop':
                break
            if startit==1:
                times.append(int(words[0]))
                pot_energies.append(float(words[1]))
            if words[0]=='Step':
                startit = 1
    infile.close()
    
    N = len(times)
    crashbefore = 0
    for i in range(N):
        pot_en   = pot_energies[i]
        thistime = times[i]
        if (pot_en!=0 and crashbefore==0):
            # Append collition times
            crashtimes.append(thistime)
            outfile.write('%i %.16f\n' % (thistime, pot_en))
            print('crashtime: %i, pot.en.:, %.5f' % (thistime, pot_en))
            crashbefore = 1
            #
            timeint       = thistime-prev_coll
            av_coll_int  += timeint # Should I bin the data instead? --> I think this is a bit murky as not all simulations will have the same number of collisions.
            prev_coll     = thistime # Or should I test when pot_en==0 again and use that as prev_coll? Nah.
            coll_int_every.append(timeint)
            coll_counter += 1
        elif pot_en==0:
            crashbefore = 0
    
    all_crashtimes.append(crashtimes)
print('%i out of %i configs skipped' % (skippedfiles, Nconfigs))
print('%i out of %i configs read' % (filescounter,Nconfigs))
outfile.close()

av_coll_int /= (coll_counter)
rms_coll_int = 0
for i in range(coll_counter):
    rms_coll_int += (coll_int_every[i]-av_coll_int)**2
rms_coll_int = np.sqrt(rms_coll_int/(coll_counter-1))

outfile_average = open(outfilename_average, 'w')
outfile_average.write('%.16f %.16f' % (av_coll_int, rms_coll_int))
outfile_average.close()
