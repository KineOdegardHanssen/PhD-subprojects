import numpy as np
import matplotlib.pyplot as plt

################################ HH #############################
# Default HH values:
ena = 50
ek = -77
el_hh = -54.3
gnabar_hh = 0.12
gkbar_hh = 0.036
gl_hh = 0.0003

### Change HH values here: ####
#ena = 45
#ek = -85
#el_hh = -40
#gnabar_hh = 0.22
#gkbar_hh = 0.020
#gl_hh = 0.001

######################### Other params ##########################
idur = 1000
somasize = 10
cms = [0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]

infolder = 'Results/IStim/Soma%i/Vary_iamp/' % somasize
hhstring = '_ena'+str(ena)+'_ek'+str(ek)+'_el'+str(el_hh)+'_gnabar'+str(gnabar_hh)+'_gkbar'+str(gkbar_hh)+'_gl'+str(gl_hh)

plotname = infolder+'somaonly_idur%i_varycm'% idur+hhstring+'_Nspikes_vs_iamp.png'

# Set names

plt.figure(figsize=(6,5))
for cm in cms:
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
    
    plt.plot(iamp, Nspikes, label=r'$C_m$=%s' % str(cm))
plt.xlabel(r'$I$ (nA)')
plt.ylabel('Number of spikes')
plt.title('Number of spikes vs input current')
plt.legend(loc='lower right')
plt.savefig(plotname)
plt.show()