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

# Bug-prone, so harcoded I:
Is     = [0,0.02,0.04,0.06,0.065,0.0655,0.066,0.0665,0.067,0.0675,0.068,0.0685,0.069,0.0695,0.07,0.0705,0.071,0.0715,0.072,0.0725,0.073,0.0735,0.074,0.0745,0.075,0.0755,0.076,0.0765,0.077,0.0775,0.078,0.0785,0.079,0.0795,0.08,0.0805,0.081,0.0815,0.082,0.0825,0.083,0.0835,0.084,0.0845,0.085, 0.09, 0.1, 0.12, 0.14, 0.16, 0.17, 0.18, 0.19, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44] 
NI     = len(Is)
cm     = 1.0 
cmfacs = [0.5,0.67,0.75,1.0,1.25,1.5,1.75,2,2.5,3]
N      = len(cmfacs)

outfolder   = 'Results/IStim/Soma%i/' % somasize
plotname    = outfolder + 'thresholds_bashh_somaprox_varyfactor'+plotstring+'.png' 
outfilename = outfolder + 'thresholds_bashh_somaprox_varyfactor'+plotstring+'.txt' 
outfile     = open(outfilename,'w')

thresholds = np.zeros(N)

for l in range(N):
    cmfac = cmfacs[l]
    changestring =''
    if varyE_bool==True:
        changestring = changestring+'_E'+str(varyE)+'_gdflt'
    else:
        changestring = changestring+'_Edefault_g'+str(varyg)
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
        if Npeaks>0:
            thresholds[l]=iamp
            outfile.write('%.2f %.5f\n' % (cmfac,iamp))
            break # break or continue?

outfile.close()

plt.figure(figsize=(6,5))
plt.plot(cmfacs,thresholds,'-o')
plt.xlabel('$C_m$ ($\mu$F/cm$^2$)')
plt.ylabel('Threshold (nA)')  
plt.title(r'Spiking threshold vs $C_m$, BAS-HH')
plt.savefig(plotname)
plt.tight_layout()
plt.show()

