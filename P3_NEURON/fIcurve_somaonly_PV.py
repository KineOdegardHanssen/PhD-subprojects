import matplotlib.pyplot as plt
import numpy as np

testmodels = [478513437]#[478513407,478513437,488462965]
cm_factors = [0.5,1,1.5]#[0.5,0.75,1,1.25,1.5,2,3,4,5] # Reduce a bit?
iamps      = np.linspace(0,0.44,23)
idur       = 1000
idelay     = 10
tstop      = idur+idelay+10. #2.
dtexp      = -8 # Soma only, so this should be fast regardless
somasize   = 10
Nmodels    = len(testmodels)
Ncm        = len(cm_factors)
Ni         = len(iamps)

# Avg and rms and stuff...
plotfolder = 'Results/IStim/Soma%i/' % somasize
plotname   =  plotfolder + 'fI_somaonlyPV_varycmfactor_%icms_idur%i.png' % (Ncm,idur)

for i in range(Ncm):
    cm_factor = cm_factors[i]
    Nspikes_avg = np.zeros(Ni)
    Nspikes_rms = np.zeros(Ni)
    for j in range(Ni):
        iamp = iamps[j]
        folder = 'Results/IStim/Soma%i/current_idur%i_iamp'%(somasize,idur)+str(iamp)+'/'
        theseNspikes = []
        for k in range(Nmodels):
            testmodel = testmodels[k]
            filename = folder+'somaonly_model%i_cmfactor'%testmodel+str(cm_factor)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_V.txt' % dtexp

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
            vthr  = vmax-0.15*deltav # If there is a peak above this value, we count it
            vprev = vthr-40 # A peak never kicks in at initiation, even if I change vthr
            Npeaks = 0
            for m in range (1,len(V)-1):
                if V[m-1]<V[m] and V[m+1]<V[m] and V[m]>vthr:
                    Npeaks+=1
            if vmax<-40: # This is not a spike # Too conservative? Should have chosen -20?
                Npeaks = 0
        
            Nspikes_avg[j]+=Npeaks
            theseNspikes.append(Npeaks)
        Nspikes_avg[j]/=Nmodels
        avg = Nspikes_avg[j]
        rms = 0
        for l in range(Nmodels):
            rms+=(avg-theseNspikes[l])**2
        if Nmodels>1:
            Nspikes_rms[j] = np.sqrt(rms/(Nmodels-1))
        else:
            Nspikes_rms[j] = 0
    if Nmodels>1:
        plt.errorbar(iamps,Nspikes_avg,yerr=Nspikes_rms,capsize=2,label=r'%s*$C_m$' % str(cm_factor))
    else:
        plt.plot(iamps,Nspikes_avg,'-o',label=r'%s*$C_m$, %i' % (str(cm_factor),testmodel))
plt.xlabel('I (nA)')
plt.ylabel('f (Spikes/s)')
plt.title(r'fI-curve, perisomatic ball-and-stick, changing $C_m$ soma and prox.dend.')
plt.legend(loc='upper left')
plt.savefig(plotname)
plt.show()