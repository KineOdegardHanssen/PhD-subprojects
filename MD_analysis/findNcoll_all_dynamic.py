import matplotlib.pyplot as plt                     # To plot
from scipy.optimize import curve_fit
from pylab import *
from scipy import ndimage                           # For Euclidean distance measurement
import data_treatment as datr
import numpy as np
import random
import math
import time
import os
import glob
import copy

## NB: This file only works with indices, so that it is easier to compare with the trajectories

long    = True
damp    = 10
spacing = 4.5
sigma_atom = 1
configs    = [1,200,500,750,999]
Nconfigs   = len(configs)
Nsteps     = 2001

basefolder = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Spacing'+str(spacing)+'/damp%i_diffseedLgv/Brush/Sigma_bead_' % damp+str(sigma_atom) + '/'
folder = basefolder + 'log-files/'
if long==True:
    folder = folder + 'long/'
outfolder = 'C:/Users/Kine/Documents/Projects_PhD/P2_PolymerMD/Planar_brush/Diffusion_bead_near_grid/Singletrajectories/Spacing' + str(spacing) + '/'



for config in configs: # Write to file inside loop?
    print('config:', config)
    
    filename_infolder = 'log.confignr%i_printevery10' % config
    if long==True:
        filename_infolder = filename_infolder + '_long'
    infilename = folder + filename_infolder
    #print('infilename:', infilename)
    
    outfilename = outfolder + 'collisiontimes_d'+str(spacing)+'_config%i.txt' % config
    outfile = open(outfilename, 'w')
    
    try:
        infile_all = open(infilename, "r")
    except:
        try:
            filename_infolder = 'log.confignr%i' % config
            infilename = folder + filename_infolder
            infile_all = open(infilename, "r")
        except:
            print('Oh, log-file! Where art thou?')
            continue # Skipping this file if it does not exist
    # Moving on, if the file
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
            crashbefore = 1
        elif pot_en==0:
            crashbefore = 0
    
    outfile.close()
