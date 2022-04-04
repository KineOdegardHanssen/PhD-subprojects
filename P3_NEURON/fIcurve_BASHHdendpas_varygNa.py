import matplotlib.pyplot as plt
import numpy as np

cm         = 1.0
iamps      = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5]
cm         = 1.0
spikedurat = -40
idur       = 1000 #100 # ms
idelay     = 100
v_init     = -65 # mV
Ra         = 100
somasize   = 10 # 15 # 
dendlen    = 1000
denddiam   = 1
nsegments  = 200 
Ni         = len(iamps)

# Default HH values:
ena = 50
ek = -77
el_hh = -54.3
gnabar_hh = 0.12
gkbar_hh = 0.036
gl_hh = 0.0003

folderstring = 'VaryNa/'
namestring = 'na'
gnabar_hh_factors = [1.0,1.2,1.5,2.0,2.5]#,3.0,5.0,10.0]
Ngs = len(gnabar_hh_factors)
varygs = np.zeros(Ngs)
for i in range(Ngs):
    varygs[i] = gnabar_hh*gnabar_hh_factors[i]

# Avg and rms and stuff...
plotfolder = 'Results/IStim/'
plotname   =  plotfolder + 'fI_BASHHdendpas_cm'+str(cm)+'_idur%i_varygNa.png' % idur

plotcolors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan','b','k','g','r','m','indigo','maroon','darkolivegreen','dodgerblue','palevioletred','darkseagreen','midnightblue','yellow','darkgoldenrod']

Nspikes = []
for i in range(Ngs):
    infolder = 'Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+folderstring
    infilename = infolder+'basHHdpas_manual_cm'+str(cm)+'_g'+namestring+str(varygs[i])+'_Nspikes_vs_Iamp.txt'
    Npeaksbefore = 0
    theseiamps   = []
    theseNspikes = []
    
    infile = open(infilename,'r')
    lines = infile.readlines()
    for line in lines:
        words = line.split()
        if len(words)>0:
            theseiamps.append(float(words[0]))
            theseNspikes.append(float(words[1]))
    infile.close()
        
    plt.plot(theseiamps,theseNspikes,color=plotcolors[i],label=r'%.1f*$\bar{g}_{\mathregular{Na}}$' % gnabar_hh_factors[i])
    Nspikes.append(theseNspikes)

plt.xlabel('I (nA)')
plt.ylabel('f (Spikes/s)')
plt.title(r'fI-curve, BAS, changing $\bar{g}_{\mathregular{Na}}$')
plt.legend(loc='lower right',ncol=2)
plt.savefig(plotname)
plt.show()