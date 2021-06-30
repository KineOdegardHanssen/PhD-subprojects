import numpy as np
import matplotlib.pyplot as plt

idur     = 1000
somasize = 10
cms      = [0.01,0.1,1.0,1.25,5.0]
Ncms     = len(cms)

currents_Tewari = np.array([0,20,40,60,80,100,120,140,160,180])*0.001

GBM14_freq = np.array([0,5,10,20,25,50,75,100,110,125])
GBM22_freq = np.array([0,7,12.5,25,40,62.5,90,107,120,140])
Sham_freq  = np.array([0,8,20,49,80,120,150,190,220,240])
ChABC_freq = np.array([0,0,0,5,20,50,75,110,150,175])


iamp_all = []
Nspikes_all = []

infolder = 'Results/IStim/Soma%i/Vary_iamp/' % somasize
plotname = infolder+'somaonly_idur%i_varycm_Nspikes_vs_iamp_withallTewari.png'% idur
plotname_norm = infolder+'somaonly_idur%i_varycm_Nspikes_vs_iamp_withallTewari_normalized.png'% idur

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
    
    iamp_all.append(iamp)
    Nspikes_all.append(Nspikes)
    
    plt.plot(iamp, Nspikes, label=r'$C_m$=%s' % str(cm))
plt.plot(currents_Tewari,GBM14_freq, label='GMB14')
plt.plot(currents_Tewari,GBM22_freq, label='GMB22')
plt.plot(currents_Tewari,Sham_freq, label='Sham')
plt.plot(currents_Tewari,ChABC_freq, label='ChABC')
plt.xlabel(r'$I$ (nA)')
plt.ylabel('Number of spikes')
plt.title('Number of spikes vs input current')
plt.legend(loc='upper left')
plt.savefig(plotname)


print('Nspikes_all:',Nspikes_all)
print('iamp_all:',iamp_all)
plt.figure(figsize=(6,5))
for i in range(Ncms):
    iamp_this = iamp_all[i]
    Nspikes_this = np.array(Nspikes_all[i])
    plt.plot(iamp_this, np.array(Nspikes_this)/max(Nspikes_this), label=r'$C_m$=%s' % str(cms[i]))
plt.plot(currents_Tewari,GBM14_freq/max(GBM14_freq), label='GMB14')
plt.plot(currents_Tewari,GBM22_freq/max(GBM22_freq), label='GMB22')
plt.plot(currents_Tewari,Sham_freq/max(Sham_freq), label='Sham')
plt.plot(currents_Tewari,ChABC_freq/max(ChABC_freq), label='ChABC')
plt.xlabel(r'$I$ (nA)')
plt.ylabel('Number of spikes/max number of spiked')
plt.title('Normalized number of spikes vs input current')
plt.legend(loc='upper left')
plt.savefig(plotname_norm)
plt.show()

# More plots, maybe...