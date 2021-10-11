import numpy as np
import matplotlib.pyplot as plt

# Write to file too?

dtexp      = -8
smallCm    = True
cms        = [0.5,1,1.75,3]#[0.5,0.67,0.75,1,1.1,1.25,1.75,2,2.5,3]#[0.5,0.67,0.75,1,1.1,1.25,1.75,2,2.5,3]#[0.5,1,1.75,3]
varymech   = 'Kd'#'NaV' #'SK'#'Im_v2'#'Kv2like'#'pas' #
varyE_bool = True
varyE      = -107 #[-130,-107,-100,-70]#[40,50,53,65,70]#[-60,-53,-30,-10,10,30,60]
varyg      = 'None'
Nvary      = len(cms)
N          = Nvary # Hacky safeguard
    
varylist = [] # Should be redundant
plotstring = '_vary'
if varyE_bool==True:
    elem = varyE
    varylist = varyE
    plotstring = plotstring + 'E'
else:
    elem = varyg
    varylist = varyg
    plotstring = plotstring + 'g'
      
if varymech=='NaV':
    folderstring = 'VaryNa/' 
    plotstring = plotstring + '_NaV'
elif varymech=='pas':
    folderstring = 'VaryPas/'
    plotstring = plotstring + '_Pas'
elif varymech=='Kd':
    folderstring = 'VaryKd/'
    plotstring = plotstring + '_Kd'
elif varymech=='Kv2like':
    folderstring = 'VaryKv2like/'
    plotstring = plotstring + '_Kv2like'
elif varymech=='Kv3_1':
    folderstring = 'VaryKv3_1/'
    plotstring = plotstring + '_Kv3_1'
elif varymech=='SK':
    folderstring = 'VarySK/'
    plotstring = plotstring + '_SK'
elif varymech=='K_T':
    folderstring = 'VaryK_T/'
    plotstring = plotstring + '_K_T'
elif varymech=='Im_v2':
    folderstring = 'VaryIm_v2/'
    plotstring = plotstring + '_Im_v2'
dtexp = -8
idur = 1000
idelay   = 100.0
afteri   = 100.0
tstart   = -500.
tstop    = idur+afteri+idelay
v_init   = -86.5
somasize = 10 # 15 # 
dendlen  = 1000
denddiam = 2

# Might need to include more Is because of the different models
Is = [0.115,0.116,0.117,0.118,0.119,0.12]#[0.21,0.211,0.212,0.213,0.214,0.215,0.22,0.225,0.23] #[0.0591,0.0592,0.0593,0.0594,0.0595,0.0596,0.0597,0.0598,0.0599,0.06,0.07,0.075,0.08,0.085,0.09,0.1,0.12]#[0.08,0.081,0.082,0.083,0.084,0.085,0.086,0.087,0.088,0.089,0.09,0.091,0.092,0.093,0.094,0.095,0.096,0.097,0.098,0.099,0.1]#[0.06,0.08,0.081,0.082,0.083,0.084,0.085,0.086,0.087,0.088,0.089,0.09,0.091,0.092,0.093,0.094,0.095,0.096,0.097,0.098,0.099,0.1]#
# Is: [0.0591,0.0592,0.0593,0.0594,0.0595,0.0596,0.0597,0.0598,0.0599,0.06,0.07,0.075,0.08,0.085,0.09,0.1,0.12]#[0.08,0.09,0.091,0.092,0.093,0.094,0.095,0.096,0.097,0.098,0.099, 0.1]#
NI = len(Is)
models  = [478513437]#[489931686]#[478513437]#[478513407]#[488462965]#[478513407,478513437,489931686,485694403,488462965] # [478513407,478513437,488462965] # 
Nmodels = len(models)

outfolder = 'Results/Soma%i/' % somasize
plotname  = outfolder + 'thresholds_baspv_somaprox_factor'+plotstring+'.png' 
outfilename  = outfolder + 'thresholds_baspv_somaprox_factor'+plotstring+'.txt' 
outfile     = open(outfilename,'w')

thresholds = np.zeros(N)

for l in range(Nvary):
    cm = cms[l]
    print('cmf:',cm, '; ', l, ' out of ', Nvary-1)
    changestring =''
    if varyE_bool==True:
        varyE = elem
        changestring = changestring+'_E'+str(varyE)+'_gdflt'
    else:
        varyg = elem
        changestring = changestring+'_Edefault_g'+str(varyg) ## Have I used this?
    outfilename_all = outfolder + 'fIcurves_baspv_somaprox_factor'+plotstring+changestring+'.txt' 
    outfile_all = open(outfilename_all,'w')
    Nspikes_avg = np.zeros(NI)
    Nspikes_rms = np.zeros(NI)
    for j in range(Nmodels):
        model = models[j]
        for i in range(NI):
            iamp = Is[i]
            if iamp==0.0:
                iamp=0
            theseNspikes = []
            currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
            infolder = 'Results/%i/IStim/Soma%i/dendlen%i/denddiam'% (model,somasize,dendlen)+str(denddiam)+'/'+ folderstring
            infolder = infolder+currentfolder
            try:
                try:
                    changestring =''
                    if varyE_bool==True:
                        varyE = elem
                        changestring = changestring+'_E'+str(varyE)+'_gdf'
                    else:
                        varyg = elem
                    infilename = infolder+'baspv_cmf'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_vinit'+str(v_init)+changestring+'_sprx_V.txt'
                    infile = open(infilename,'r')
                except:
                    infilename = infolder+'baspv_cmf'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+changestring+'_sprx_V.txt'
                    infile = open(infilename,'r')
                    
            except:
                try:
                    infilename = infolder+'baspv_cmf'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+changestring+'_sprx_V.txt' # Scary
                    infile = open(infilename,'r')
                except:
                    changestring =''
                    if varyE_bool==True:
                        varyE = elem
                        changestring = changestring+'_E'+str(varyE)+'_gdflt'
                    else:
                        varyg = elem
                    infile = open(infilename,'r')
                
            outfilename = infolder+'baspv_cmf'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+changestring+'_sprx_If.txt'
            
            t = []
            V = []
            
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
            vthr  = vmax-0.15*deltav # If there is a peak above this value, we count it # OLD
            vprev = vthr-40 # A peak never kicks in at initiation, even if I change vthr
            Npeaks = 0
            for m in range (1,len(V)-1):
                if V[m-1]<V[m] and V[m+1]<V[m] and V[m]>vthr:
                    Npeaks+=1
            if vmax<-40: # This is not a spike # Too conservative? Should have chosen -20?
                Npeaks = 0
            print('iamp:',iamp, '; ', i, ' out of ', NI-1, '; Npeaks:', Npeaks)
            if Npeaks>0:
                thresholds[l]=iamp
                outfile.write('%.2f %.5f\n' % (cm,iamp))
                break # break or continue?

plt.figure(figsize=(6,5))            
plt.plot(cms,thresholds,'-o')
plt.xlabel('$C_m$ ($\mu$F/cm$^2$)')
plt.ylabel('Threshold (nA)')  
plt.title(r'Spiking threshold vs $C_m$, BAS-PV')
plt.savefig(plotname)
plt.show()