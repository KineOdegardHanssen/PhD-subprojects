import matplotlib.pyplot as plt                     # To plot
from pylab import *
import numpy as np
import math

Nchains = 9
Nbeads  = 100
hsub    = 1
hcyl    = 100
r       = 0.5
sigma   = 1

# Find Lz
d = np.array([1,1.25,1.5,2,3,4,5,6,7,10,15,25,50,75,100])
N = len(d)

def give_phi_nooverlap(d,sigma):
    return 1 - np.pi*sigma**2/d**2

def give_phi_partialoverlap(d,sigma):
    return 18*(sigma**2*math.acos(d/(2*sigma))-d/4.*np.sqrt(4*sigma**2-d**2))


phi = np.zeros(N)
# Total overlap:
phi[0] = 0
phi[1] = 0
# Some overlap (but ensure that three or more circles do not occupy the same area):
phi[2] = give_phi_partialoverlap(d[2], sigma)
# No overlap:
for i in range(3,N):
    phi[i] = give_phi_nooverlap(d[i],sigma)

### Plotting:

filename = 'd_vs_phi_fromflat.txt'
file = open(filename, 'w')
file.write('d phi\n')
for i in range(len(d)):
    file.write('%.2f %.16f\n' % (d[i], phi[i]))
file.close()
