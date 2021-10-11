import matplotlib.pyplot as plt
import numpy as np

cms        = [1.0,1.25]#[0.67,1.0]#,1.5]#[0.5,0.75,1,1.25,1.5,2,3,4,5] # Reduce a bit?
iamps      = np.linspace(0,0.44,23)
idur       = 1000
idelay     = 10
tstop      = idur+idelay+10. #2.
dtexp      = -8 # Soma only, so this should be fast regardless
somasize   = 10
Ra         = 100.
v_init     = -70 
Ncm        = len(cms)
Ni         = len(iamps)

# Default HH values:
ena = 50
ek = -77
el_hh = -54.3
gnabar_hh = 0.12
gkbar_hh = 0.036
gl_hh = 0.0003

hhstring = '_ena'+str(ena)+'_ek'+str(ek)+'_el'+str(el_hh)+'_gnabar'+str(gnabar_hh)+'_gkbar'+str(gkbar_hh)+'_gl'+str(gl_hh)

# Avg and rms and stuff...
plotfolder = 'Results/IStim/Soma%i/' % somasize
plotname   =  plotfolder + 'fI_somaonlyHH_varycmfactor_%icms_idur%i' % (Ncm,idur)+hhstring+'.png'

Nspikes = []
fs_at_180pA = []
fs_at_440pA = []
for i in range(Ncm):
    cm = cms[i]
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
        for m in range (1,len(V)-1):
            if V[m-1]<V[m] and V[m+1]<V[m] and V[m]>vthr:
                Npeaks+=1
        if vmax<-40: # This is not a spike # Too conservative? Should have chosen -20?
            Npeaks = 0
    
        theseNspikes.append(Npeaks)
        
        if iamp==0.18:
            fs_at_180pA.append(Npeaks)
        elif iamp==0.44:
            fs_at_440pA.append(Npeaks)
        
    plt.plot(iamps,theseNspikes,'-o',label=r'%s*$C_m$' % (str(cm)))
    Nspikes.append(theseNspikes)

delta_f180 = fs_at_180pA[0]-fs_at_180pA[-1]
delta_f440 = fs_at_440pA[0]-fs_at_440pA[-1]
reldiff_f180 = delta_f180/fs_at_180pA[0]
reldiff_f440 = delta_f440/fs_at_440pA[0]

print('delta_f180:',delta_f180)
print('reldiff_f180:',reldiff_f180)
print('delta_f440:',delta_f440)
print('reldiff_f440:',reldiff_f440)


plt.xlabel('I (nA)')
plt.ylabel('f (Spikes/s)')
plt.title(r'fI-curve, HH one comp, changing $C_m$')
plt.legend(loc='upper left')
plt.savefig(plotname)
plt.show()