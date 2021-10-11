import numpy as np
import matplotlib.pyplot as plt

changeconductance = 'k' # 'n' # 'l' # 

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
    conductances = [gnabar_hh,0.14,0.16,0.18,0.2,0.22,gnabar_hh*2,gnabar_hh*4,gnabar_hh*8]
elif changeconductance=='k':
    conductances = [gkbar_hh/4,gkbar_hh/3,gkbar_hh/2,0.02,0.025,0.03,gkbar_hh]
elif changeconductance=='l':
    conductances = [0,0.000003,0.00003,gl_hh,0.0004,0.0005,0.0006,0.0008,0.001]

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
    #plt.plot(iamp, Nspikes, '-o', label=r'$g$=%s' % str(conductance)) # For testing

if changeconductance=='n':
    string = '$\bar{g}_{Na}$'
elif changeconductance=='k':
    string = '$\bar{g}_{K}$'
elif changeconductance=='l':
    string ='$g_l$'
plt.xlabel(r'$I$ (nA)')
plt.ylabel('Number of spikes (frequency)')
if changeconductance=='n':
    plt.title(r'Number of spikes vs $I$ for $\bar{g}_{Na}$ ')
elif changeconductance=='k':
    plt.title(r'Number of spikes vs $I$ for $\bar{g}_{K}$ ')
elif changeconductance=='l':
    plt.title(r'Number of spikes vs $I$ for $g_l$ ')
plt.legend(loc='lower right')
plt.savefig(plotname)
plt.show()