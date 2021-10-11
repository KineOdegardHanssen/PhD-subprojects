import numpy as np
import matplotlib.pyplot as plt

# Write to file too?

dtexp = -8
smallCm = True
#if smallCm==True:
#    cms = [0.5,0.75,1.0,1.25,1.5,2.0] #[0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0]
#else:
#    cms = [0.5,1,2,3,4,4.5,5,5.5,6,7,8,9,10]
cms = [0.5,1,1.5]
idur = 1000
idelay   = 100.0
afteri   = 100.0
tstart   = -500.
tstop    = idur+afteri+idelay
v_init   = -86.5
somasize = 10 # 15 # 
dendlen  = 1000
denddiam = 2

NI   = 21 # 10 # 
Imax = 0.18
if NI==21:
    Imax = 0.4
Is   = np.linspace(0,Imax,NI)
models  = [478513437,488462965] #  [478513407,478513437,488462965] # 
Nmodels = len(models)
Ncms    = len(cms)

outfolder = 'Results/Soma%i/' % somasize
plotname  = outfolder + 'fIcurves_basps_somaprox_changelinearly.png' 

plt.figure(figsize=(6,5))
for l in range(Ncms):
    cm = cms[l]
    Nspikes_avg = np.zeros(NI)
    Nspikes_rms = np.zeros(NI)
    for i in range(NI):
        iamp = Is[i]
        if iamp==0.0:
            iamp=0
        if cm==1:
            cm = 1.0
        theseNspikes = []
        currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
        for j in range(Nmodels):
            model = models[j]
            infolder = 'Results/%i/IStim/Soma%i/dendlen%i/denddiam'% (model,somasize,dendlen)+str(denddiam)+'/'
            infolder = infolder+currentfolder
            infilename = infolder+'basps_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_vinit'+str(v_init)+'_somaprox_V.txt' 
            
            t = []
            V = []
            
            infile = open(infilename,'r')
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
            print(Npeaks, ' peaks for model ', model, ', current ', iamp)
            print('vmax:', vmax)
            
            Nspikes_avg[i]+=Npeaks
            theseNspikes.append(Npeaks)
        Nspikes_avg[i]/=Nmodels
        avg = Nspikes_avg[i]
        rms = 0
        for k in range(Nmodels):
            rms+=(avg-theseNspikes[k])**2
        #Nspikes_rms[i] = np.sqrt(rms/Nmodels)
        Nspikes_rms[i] = np.sqrt(rms/(Nmodels-1))
    plt.errorbar(Is,Nspikes_avg,yerr=Nspikes_rms,capsize=2,label=r'$C_m$=%s' % str(cm)) # Flexible in case we want more decimals
plt.xlabel('I (nA)')
plt.ylabel('f (Spikes/s)')
plt.title(r'fI-curve, perisomatic ball-and-stick, changing $C_m$ soma and prox.dend.')
plt.legend(loc='upper left')
plt.savefig(plotname)
plt.show()