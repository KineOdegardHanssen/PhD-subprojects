import matplotlib.pyplot as plt
import numpy as np

cm         = 1.0
iamps      = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5]
idur       = 1000
idelay     = 10
tstop      = idur+idelay+10. #2.
dtexp      = -8 # Soma only, so this should be fast regardless
somasize   = 10
Ra         = 100.
v_init     = -70 
Ni         = len(iamps)

# Default HH values:
ena = 50
ek = -77
el_hh = -54.3
gnabar_hh = 0.12
gkbar_hh = 0.036
gl_hh = 0.0003

gnabar_hh_factors = [0.7,0.75,0.8,0.9,1.0,1.1,1.2,1.5]#[0.7,0.8,0.9,1.0,1.1,1.2,1.5,2,3,4,5,10]
Ngs = len(gnabar_hh_factors)
gnabars_hh = np.zeros(Ngs)
for i in range(Ngs):
    gnabars_hh[i] = gnabar_hh*gnabar_hh_factors[i]

# Avg and rms and stuff...
plotfolder = 'Results/IStim/'
plotname   =  plotfolder + 'fI_somaonlyHH_varysize_cm'+str(cm)+'_idur%i_varygNa.png' % idur

plotcolors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan','b','k','g','r','m','indigo','maroon','darkolivegreen','dodgerblue','palevioletred','darkseagreen','midnightblue','yellow','darkgoldenrod']

Nspikes = []
for i in range(Ngs):
    gnabar_hh = gnabars_hh[i]
    hhstring = '_ena'+str(ena)+'_ek'+str(ek)+'_el'+str(el_hh)+'_gnabar'+str(gnabar_hh)+'_gkbar'+str(gkbar_hh)+'_gl'+str(gl_hh)
    Npeaksbefore = 0
    theseiamps   = []
    theseNspikes = []
    for j in range(Ni):
        iamp = iamps[j]
        folder = 'Results/IStim/Soma%i/current_idur%i_iamp'%(somasize,idur)+str(iamp)+'/'
        filename = folder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+hhstring+'_Ra%i_vinit' %Ra+str(v_init)+'_V.txt' 

        t = []
        V = []
        
        infile = open(filename,'r')
        lines = infile.readlines()
        for line in lines:
            words = line.split()
            if len(words)>0:
                t.append(float(words[0]))
                V.append(float(words[1]))
        infile.close()
        
        vmax = max(V) 
        vmin = min(V) 
        deltav = vmax-vmin
        vthr  = -40 #vmax-0.15*deltav # If there is a peak above this value, we count it
        vprev = vthr-40 # A peak never kicks in at initiation, even if I change vthr
        Npeaks = 0
        lastpeaktime = 8000
        peaktimes = []
        for m in range (1,len(V)-1):
            if V[m-1]<V[m] and V[m+1]<V[m] and V[m]>vthr:
                Npeaks+=1
                lastpeaktime = t[m]
                peaktimes.append(lastpeaktime)
        if vmax<-40: # This is not a spike # Too conservative? Should have chosen -20?
            Npeaks = 0
        if lastpeaktime<idelay+(idur/2.):
            Npeaks = 0
        middlestim = 0
        for k in range(Npeaks):
            if peaktimes[k]<idelay+0.8*idur and peaktimes[k]>idelay+0.2*idur:
                middlestim+=1
        if middlestim==0: # No peak in the middle of the stimulation interval
            Npeaks = 0
        if Npeaksbefore!=0 and Npeaks==0:
            break
        Npeaksbefore = Npeaks
        theseiamps.append(iamp)
        theseNspikes.append(Npeaks)
        outfilename = folder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+hhstring+'_Ra%i_vinit' %Ra+str(v_init)+'_Nspikes_vs_I.txt'
        outfile = open(outfilename,'w')
        outfile.write('%i' % Npeaks)
        outfile.close()
        
    plt.plot(theseiamps,theseNspikes,color=plotcolors[i],label=r'%s*$\bar{g}_{\mathregular{Na}}$' % str(gnabar_hh_factors[i]))
    Nspikes.append(theseNspikes)

plt.xlabel('I (nA)')
plt.ylabel('f (Spikes/s)')
plt.title(r'fI-curve, HH one comp, changing $C_m$')
plt.legend(loc='lower right',ncol=2)
plt.savefig(plotname)
plt.show()