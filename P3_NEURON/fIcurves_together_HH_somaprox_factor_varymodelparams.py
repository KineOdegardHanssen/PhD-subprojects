import numpy as np
import matplotlib.pyplot as plt

# Write to file too?

dtexp = -8
smallCm = True
cm = 1
varymech = 'Na' #'K'#'leak'
varyE_bool = True
varyE = [-10,20,50,60]#[-130,-107,-100,-70]#[40,50,53,65,70]#[-60,-53,-30,-10,10,30,60]
varyg = 'None'
    
varylist = [] # Should be redundant
plotstring = '_vary'
if varyE_bool==True:
    varylist = varyE
    plotstring = plotstring + 'E'
else:
    varylist = varyg
    plotstring = plotstring + 'g'
Nvary    = len(varylist)
      
if varymech=='Na':
    folderstring = 'VaryNa/' 
    plotstring = plotstring + '_Na'
elif varymech=='pas':
    folderstring = 'VaryLeak/'
    plotstring = plotstring + '_leak'
elif varymech=='K':
    folderstring = 'VaryK/'
    plotstring = plotstring + '_K'

dtexp = -8
idur = 1000
idelay   = 100.0
afteri   = 100.0
tstart   = -500.
tstop    = idur+afteri+idelay
v_init   = -65
Ra       = 100
somasize = 10 # 15 # 
dendlen  = 1000
denddiam = 2

''' # Bug-prone, so will harcode instead :(
NI   = 23 # 10 # 
Imax = 0.18
if NI==23:
    Imax = 0.44
Is   = np.linspace(0,Imax,NI)
'''
Is = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44] 
NI = len(Is)
cm        = 1.0 
cmfac     = 1.25

outfolder = 'Results/IStim/Soma%i/' % somasize
plotname  = outfolder + 'fIcurves_baspv_somaprox_factor'+plotstring+'.png' 


plt.figure(figsize=(6,5))
for l in range(Nvary):
    elem = varylist[l]
    changestring =''
    if varyE_bool==True:
        varyE = elem
        changestring = changestring+'_E'+str(varyE)+'_gdflt'
    else:
        varyg = elem
        changestring = changestring+'_Edefault_g'+str(varyg)
    outfilename_all = outfolder + 'fIcurves_basHH_somaprox_factor'+plotstring+changestring+'.txt' 
    outfile_all = open(outfilename_all,'w')
    Nspikes     = np.zeros(NI)
    for i in range(NI):
        iamp = Is[i]
        if iamp==0.0:
            iamp=0
        theseNspikes = []
        currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
        infolder = 'Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+ folderstring
        infolder = infolder+currentfolder
        infilename = infolder+'basHH_cmf'+str(cmfac)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+changestring+'_V_sprxf.txt' 
        outfilename = infolder+'basHH_cmf'+str(cmfac)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+changestring+'_V_sprxf_If.txt' 
        
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
        vthr  = -40#vmax-0.15*deltav # If there is a peak above this value, we count it # OLD
        vprev = vthr-40 # A peak never kicks in at initiation, even if I change vthr
        Npeaks = 0
        for m in range (1,len(V)-1):
            if V[m-1]<V[m] and V[m+1]<V[m] and V[m]>vthr:
                Npeaks+=1
        if vmax<-40: # This is not a spike # Too conservative? Should have chosen -20?
            Npeaks = 0
        print(Npeaks, ' peaks, current ', iamp)
        print('vmax:', vmax)
        
        outfile = open(outfilename,'w')
        outfile.write('%.2f %i' % (iamp,Npeaks))
        outfile.close()
        
        outfile_all.write('%.2f %i\n' % (iamp,Npeaks))
        Nspikes[i]=Npeaks
    if varyE_bool==True:
        plt.plot(Is,Nspikes,label=r'E=%s' % str(elem)) # Flexible in case we want more decimals
    elif varyE_bool==False: # Might want to add plotting for different Cm's
        plt.plot(Is,Nspikes,label=r'E=%s' % str(elem)) # Flexible in case we want more decimals
    outfile_all.close()
plt.xlabel('I (nA)')
plt.ylabel('f (Spikes/s)')      
if varymech=='NaV':
    if varyE_bool==True:
        plt.title(r'fI-curve, perisomatic ball-and-stick, changing $E_{Na}$ soma and prox.dend.')
    else:
        plt.title(r'fI-curve, perisomatic ball-and-stick, changing $\bar{g}_{Na}$ soma and prox.dend.')
plt.legend(loc='upper left')
plt.savefig(plotname)
plt.show()