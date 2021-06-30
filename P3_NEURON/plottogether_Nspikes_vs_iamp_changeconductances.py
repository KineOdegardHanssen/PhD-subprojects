import numpy as np
import matplotlib.pyplot as plt

changeconductance = 'l' # 'n' # 'k' # 

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

if changeconductance=='n': 
    conductances = [gnabar_hh,gnabar_hh*2,gnabar_hh*4,gnabar_hh*8]
elif changeconductance=='k':
    conductances = [gkbar_hh/4,gkbar_hh/3,gkbar_hh/2,gkbar_hh]
elif changeconductance=='l':
    conductances = [0,0.000003,0.00003,gl_hh]

infolder = 'Results/IStim/Soma%i/Vary_iamp/' % somasize
## If-test for name here too

if changeconductance=='n':
    plotname = infolder+'somaonly_idur%i_cm'% idur+str(cm)+'_varyhhparam_gnabar_Nspikes_vs_iamp.png'
elif changeconductance=='k':
    plotname = infolder+'somaonly_idur%i_cm'% idur+str(cm)+'_varyhhparam_gkbar_Nspikes_vs_iamp.png'
elif changeconductance=='l':
    plotname = infolder+'somaonly_idur%i_cm'% idur+str(cm)+'_varyhhparam_gl_Nspikes_vs_iamp.png'

plt.figure(figsize=(6,5))
for conductance in conductances:
    if changeconductance=='n':
        gnabar_hh = conductance
    elif changeconductance=='k':
        gkbar_hh = conductance
    elif changeconductance=='l':
        gl_hh = conductance
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
    
    plt.plot(iamp, Nspikes, label=r'$g$=%s' % str(conductance))

if changeconductance=='n':
    plt.xlabel(r'$\bar{g}_{Na}$')
elif changeconductance=='k':
    plt.xlabel(r'$\bar{g}_{K}$')
elif changeconductance=='l':
    plt.xlabel(r'$g_l$')
plt.ylabel('Number of spikes (frequency)')
plt.title('Number of spikes vs g')
plt.legend(loc='lower right')
plt.savefig(plotname)
plt.show()