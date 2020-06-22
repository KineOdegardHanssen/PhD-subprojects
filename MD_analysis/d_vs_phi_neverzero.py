import matplotlib.pyplot as plt                     # To plot
from pylab import *
import numpy as np

Nchains = 9
Nbeads  = 100
hsub    = 1
hcyl    = 100
r       = 1 #0.5
Lz      = 300

def give_phiblock_wsubstr(d):
    porosity = 1 - (pi*r**2*hcyl+d**2*hsub)/(d**2*Lz)
    Np = len(porosity)
    for i in range(Np):
        if porosity[i]<0:
            porosity[i] = 0
    return porosity

def give_phibead_wsubstr(d):
    porosity = 1 - (4./3*pi*r**3*(Nbeads+0.5*d**2))/(d**2*Lz)
    Np = len(porosity)
    for i in range(Np):
        if porosity[i]<0:
            porosity[i] = 0
    return porosity

def give_phiblock_chainsonly(d):
    porosity = 1 - (pi*r**2*hcyl)/(d**2*Lz)
    Np = len(porosity)
    for i in range(Np):
        if porosity[i]<0:
            porosity[i] = 0
    return porosity

def give_phibead_chainsonly(d):
    porosity = 1 - (4./3*pi*r**3*Nbeads)/(d**2*Lz)
    Np = len(porosity)
    for i in range(Np):
        if porosity[i]<0:
            porosity[i] = 0
    return porosity

d = np.linspace(1,100,201)

phi_block = give_phiblock_wsubstr(d)
phi_bead  = give_phibead_wsubstr(d)

plotname = 'd_vs_phi_r%.1f_neverzero.png' % r

plt.figure()
plt.plot(d,phi_block, label='Block model')
plt.plot(d,phi_bead, label='Bead model')
plt.xlabel(r'Spacing $d$')
plt.ylabel(r'Porosity $\phi$')
plt.title('Spacing $d$ vs porosity $d$')
plt.tight_layout()
plt.legend(loc='lower right')
#plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.savefig(plotname)


## Short
d = np.linspace(1,20,201)

phi_block = give_phiblock_wsubstr(d)
phi_bead  = give_phibead_wsubstr(d)

plotname = 'd_vs_phi_short_r%.1f_neverzero.png' % r

plt.figure()
plt.plot(d,phi_block, label='Block model')
plt.plot(d,phi_bead, label='Bead model')
plt.xlabel(r'Spacing $d$')
plt.ylabel(r'Porosity $\phi$')
plt.title('Spacing $d$ vs porosity $d$')
plt.tight_layout()
plt.legend(loc='lower right')
#plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.savefig(plotname)

print('phi_block(d=1):', phi_block[0])
print('phi_bead(d=1):', phi_bead[0])

################# Without substrate ####################
d = np.linspace(1,100,201)

phi_block = give_phiblock_chainsonly(d)
phi_bead  = give_phibead_chainsonly(d)

plotname = 'd_vs_phi_wosubstrate_r%.1f_neverzero.png'  % r

plt.figure()
plt.plot(d,phi_block, label='Block model')
plt.plot(d,phi_bead, label='Bead model')
plt.xlabel(r'Spacing $d$')
plt.ylabel(r'Porosity $\phi$')
plt.title('Spacing $d$ vs porosity $d$')
plt.tight_layout()
plt.legend(loc='lower right')
#plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.savefig(plotname)

d = np.linspace(1,20,201)

phi_block = give_phiblock_chainsonly(d)
phi_bead  = give_phibead_chainsonly(d)

plotname = 'd_vs_phi_short_wosubstrate_r%.1f_neverzero.png' % r 

plt.figure()
plt.plot(d,phi_block, label='Block model')
plt.plot(d,phi_bead, label='Bead model')
plt.xlabel(r'Spacing $d$')
plt.ylabel(r'Porosity $\phi$')
plt.title('Spacing $d$ vs porosity $d$')
plt.tight_layout()
plt.legend(loc='lower right')
#plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.savefig(plotname)

print('phi_block(d=1):', phi_block[0])
print('phi_bead(d=1):', phi_bead[0])


### Plotting:
d = np.array([1,1.25,1.5,2,3,4,5,6,7,10,15,25,50,75,100])
phi_block_ws = give_phiblock_wsubstr(d)
phi_bead_ws  = give_phibead_wsubstr(d)
phi_block_co = give_phiblock_chainsonly(d)
phi_bead_co  = give_phibead_chainsonly(d)
filename_ws = 'd_vs_phi_wsubstrate_r%.1f_neverzero.txt' % r
file_ws = open(filename_ws, 'w')
file_ws.write('d phi(blockmodel) phi(beadmodel)\n')
for i in range(len(d)):
    file_ws.write('%.2f %.16f %.16f\n' % (d[i], phi_block_ws[i], phi_bead_ws[i]))
file_ws.close()

filename_wos = 'd_vs_phi_withoutsubstrate_r%.1f_neverzero.txt' % r
file_wos = open(filename_wos, 'w')
file_wos.write('d phi(blockmodel) phi(beadmodel)\n')
for i in range(len(d)):
    file_wos.write('%.2f %.16f %.16f\n' % (d[i], phi_block_co[i], phi_bead_co[i]))
file_wos.close()