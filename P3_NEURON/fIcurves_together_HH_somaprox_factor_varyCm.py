import numpy as np
import matplotlib.pyplot as plt

read_off = True # False # 

dtexp = -8
varymech = 'Na' #'K'#'leak'
varyE_bool = True
varyE = 50 #[-10,20,50,60]#[-130,-107,-100,-70]#[40,50,53,65,70]#[-60,-53,-30,-10,10,30,60]
varyg = 'None'
    
varylist = [] # Should be redundant
plotstring = '_vary'
if varyE_bool==True:
    varylist = varyE
    plotstring = plotstring + 'E'+str(varyE)
else:
    varylist = varyg
    plotstring = plotstring + 'g'+str(varyg)
      
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

# Bug-prone, so will harcoded I:
Is     = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44] 
NI     = len(Is)
cm     = 1.0 
if read_off==True:
    cmfacs = [0.67,1.0,1.5]#[1.0,1.25,1.5]
else:
    cmfacs = [0.1,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2,3]
N      = len(cmfacs)

outfolder = 'Results/IStim/Soma%i/' % somasize
plotname  = outfolder + 'fIcurves_bashh_somaprox_varyfactor'+plotstring+'.png' 

f_end = np.zeros(N)
f_018 = np.zeros(N)

plt.figure(figsize=(6,5))
for l in range(N):
    cmfac = cmfacs[l]
    changestring =''
    if varyE_bool==True:
        changestring = changestring+'_E'+str(varyE)+'_gdflt'
    else:
        changestring = changestring+'_Edefault_g'+str(varyg)
    outfilename_all = outfolder + 'fIcurves_basHH_somaprox_varyCmf'+plotstring+changestring+'.txt' 
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
        vthr  = -40 #vmax-0.15*deltav # If there is a peak above this value, we count it # OLD
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
        
        if read_off==True:
            if iamp==0.18:
                f_018[l] = Npeaks
        
    plt.plot(Is,Nspikes,label=r'$C_{m,factor}$=%s' % str(cmfac))
    outfile_all.close()
    
    if read_off==True:
        f_end[l] = Nspikes[-1]

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

if read_off==True:
    #deltaf25 = (f_end[1]-f_end[0])/f_end[0]
    #deltaf50 = (f_end[2]-f_end[0])/f_end[0]
    deltaf25 = (f_end[0]-f_end[1])/f_end[1]
    deltaf50 = (f_end[0]-f_end[2])/f_end[2]
    print('delta f, 25% diff in Cm:', deltaf25)
    print('delta f, 50% diff in Cm:', deltaf50)
    diff1 = f_end[0]-f_end[1]
    diff2 = f_end[0]-f_end[2]
    print('diff 1, no. of peaks:', diff1)
    print('diff 2, no. of peaks:', diff2)
    #
    deltaf25 = (f_018[0]-f_018[1])/f_018[1]
    deltaf50 = (f_018[0]-f_018[2])/f_018[2]
    print('delta f, I=0.18nA, 25% diff in Cm:', deltaf25)
    print('delta f, I=0.18nA, 50% diff in Cm:', deltaf50)
    diff1 = f_018[0]-f_018[1]
    diff2 = f_018[0]-f_018[2]
    print('diff 1, no. of peaks, I=0.18nA:', diff1)
    print('diff 2, no. of peaks, I=0.18nA:', diff2)