import matplotlib.pyplot as plt                     # To plot
from pylab import *
import numpy as np

Nchains = 9
Nbeads  = 100
hsub    = 1
hcyl    = 100
r       = 0.5
Lz      = 300

d = np.linspace(1,100,201)

phi_block = 1 - (pi*r**2*hcyl+d**2*hsub)/(d**2*Lz)
phi_bead  = 1 - (4./3*pi*r**3*(Nbeads+0.5*d**2))/(d**2*Lz)

plotname = 'd_vs_phi.png'
filename_ws = 'd_vs_phi_wsubstrate.txt'
file_ws = open(filename_ws, 'w')
file_ws.write('d phi(blockmodel) phi(beadmodel)\n')
for i in range(len(d)):
    file_ws.write('%.2f %.16f %.16f\n' % (d[i], phi_block[i], phi_bead[i]))
file_ws.close()

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

phi_block = 1 - (pi*r**2*hcyl+d**2*hsub)/(d**2*Lz)
phi_bead  = 1 - (4./3*pi*r**3*(Nbeads+0.5*d**2))/(d**2*Lz)

plotname = 'd_vs_phi_short.png'

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

phi_block = 1 - (pi*r**2*hcyl)/(d**2*Lz)
phi_bead  = 1 - (4./3*pi*r**3*Nbeads)/(d**2*Lz)

plotname = 'd_vs_phi_wosubstrate.png'
filename_wos = 'd_vs_phi_withoutsubstrate.txt'
file_wos = open(filename_ws, 'w')
file_wos.write('d phi(blockmodel) phi(beadmodel)\n')
for i in range(len(d)):
    file_wos.write('%.2f %.16f %.16f\n' % (d[i], phi_block[i], phi_bead[i]))
file_wos.close()

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

phi_block = 1 - (pi*r**2*hcyl)/(d**2*Lz)
phi_bead  = 1 - (4./3*pi*r**3*Nbeads)/(d**2*Lz)

plotname = 'd_vs_phi_short_wosubstrate.png'

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

