import numpy as np
import matplotlib.pyplot as plt

changepotential = 'n' # 'l' # 'k' # 
idur = 1000
somasize = 10
cm = 1.0

# Default HH values:
ena = 50
ek = -77
el_hh = -54.3
gnabar_hh = 0.12
gkbar_hh = 0.036
gl_hh = 0.0003

if changepotential=='n': 
    potentials = [40,45,48,49,49.3,49.5,50,60]
elif changepotential=='k':
    potentials = [-70,-77,-85,-90,-100]
elif changepotential=='l':
    potentials = [-70,-60,-54.3,-40,-30,-20,0,20]

infolder = 'Results/IStim/Soma%i/Vary_iamp/' % somasize
## If-test for name here too

if changepotential=='n':
    plotname = infolder+'somaonly_idur%i_cm'% idur+str(cm)+'_varyhhparam_Ena_Nspikes_vs_iamp.png'
elif changepotential=='k':
    plotname = infolder+'somaonly_idur%i_cm'% idur+str(cm)+'_varyhhparam_Ek_Nspikes_vs_iamp.png'
elif changepotential=='l':
    plotname = infolder+'somaonly_idur%i_cm'% idur+str(cm)+'_varyhhparam_El_Nspikes_vs_iamp.png'

plt.figure(figsize=(6,5))
for potential in potentials:
    if changepotential=='n':
        ena = potential
    elif changepotential=='k':
        ek = potential
    elif changepotential=='l':
        el_hh = potential
    hhstring = '_ena'+str(ena)+'_ek'+str(ek)+'_el'+str(el_hh)+'_gnabar'+str(gnabar_hh)+'_gkbar'+str(gkbar_hh)+'_gl'+str(gl_hh)
    infilename = infolder+'somaonly_idur%i_cm'% idur+str(cm)+hhstring+'_Nspikes_vs_iamp.txt'
    
    infile = open(infilename,'r')
    lines = infile.readlines()
    
    iamp = []
    Nspikes = []
    
    for line in lines:
        words = line.split()
        if len(words)>0:
            iamp.append(float(words[0]))
            Nspikes.append(float(words[1]))
    infile.close()
    
    plt.plot(iamp, Nspikes, label=r'$E$=%s' % str(potential))
    #plt.plot(iamp, Nspikes, '-o', label=r'$E$=%s' % str(potential))

plt.xlabel(r'$I$ (nA)')
plt.ylabel('Number of spikes (frequency)')
if changepotential=='n':
    plt.title(r'Number of spikes vs $E_{Na}$')
elif changepotential=='k':
    plt.title(r'Number of spikes vs $E_K$')
elif changepotential=='l':
    plt.title(r'Number of spikes vs $E_l$')
plt.legend(loc='lower right')
plt.savefig(plotname)
plt.show()