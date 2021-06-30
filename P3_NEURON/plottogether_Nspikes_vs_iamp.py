import numpy as np
import matplotlib.pyplot as plt

idur = 1000
somasize = 15
cms = [1.0,1.25]

infolder = 'Results/IStim/Soma%i/Vary_iamp/' % somasize
plotname = infolder+'somaonly_idur%i_varycm_Nspikes_vs_iamp.png'% idur

plt.figure(figsize=(6,5))
for cm in cms:
    infilename = infolder+'somaonly_idur%i_cm'% idur+str(cm) +'_Nspikes_vs_iamp.txt'
    
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