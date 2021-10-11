import numpy as np
import matplotlib.pyplot as plt

# Write to file too?

dtexp = -8
smallCm = True
cm = 1
varymech = 'Kd'#'NaV' #'SK'#'Im_v2'#'Kv2like'#'pas' #
varyE_bool = True
varyE = [-130,-107,-100,-70]#[40,50,53,65,70]#[-60,-53,-30,-10,10,30,60]
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
#testmodel = 488462965 # 478513437#478513407# 485694403 # 489931686 #  ### Weirdo:480633479#  
idur = 1000
idelay   = 100.0
afteri   = 100.0
tstart   = -500.
tstop    = idur+afteri+idelay
v_init   = -86.5
somasize = 10 # 15 # 
dendlen  = 1000
denddiam = 2

NI   = 23 # 10 # 
Imax = 0.18
if NI==23:
    Imax = 0.44
Is   = np.linspace(0,Imax,NI)
models  = [478513407] # [478513437] # ## [489931686] # [485694403] # [488462965] # [478513407,478513437,488462965] # 
Nmodels = len(models)

outfolder = 'Results/Soma%i/' % somasize
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
        changestring = changestring+'_Edefault_g'+str(varyg) ## Have I used this?
    outfilename_all = outfolder + 'fIcurves_baspv_somaprox_factor'+plotstring+changestring+'.txt' 
    outfile_all = open(outfilename_all,'w')
    Nspikes_avg = np.zeros(NI)
    Nspikes_rms = np.zeros(NI)
    for i in range(NI):
        iamp = Is[i]
        if iamp==0.0:
            iamp=0
        theseNspikes = []
        currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
        for j in range(Nmodels):
            model = models[j]
            infolder = 'Results/%i/IStim/Soma%i/dendlen%i/denddiam'% (model,somasize,dendlen)+str(denddiam)+'/'+ folderstring
            infolder = infolder+currentfolder
            try:
                infilename = infolder+'baspv_cmf'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_vinit'+str(v_init)+changestring+'_sprx_V.txt'
                infile = open(infilename,'r')
            except:
                if varyE_bool==True:
                    varyE = elem
                    changestring = changestring+'_E'+str(varyE)+'_gdf'
                else:
                    varyg = elem
                    changestring = changestring+'_Edf_g'+str(varyg)
                    infilename = infolder+'baspv_cmf'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_vinit'+str(v_init)+changestring+'_sprx_V.txt'
                    infile = open(infilename,'r')
                
            outfilename = infolder+'baspv_cmf'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_vinit'+str(v_init)+changestring+'_sprx_If.txt'
            
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
            print(Npeaks, ' peaks for model ', model, ', current ', iamp)
            print('vmax:', vmax)
            
            outfile = open(outfilename,'w')
            outfile.write('%.2f %i' % (iamp,Npeaks))
            outfile.close()
            
            Nspikes_avg[i]+=Npeaks
            theseNspikes.append(Npeaks)
        Nspikes_avg[i]/=Nmodels
        avg = Nspikes_avg[i]
        rms = 0
        for k in range(Nmodels):
            rms+=(avg-theseNspikes[k])**2
        #Nspikes_rms[i] = np.sqrt(rms/Nmodels)
        Nspikes_rms[i] = np.sqrt(rms/(Nmodels-1))
        outfile_all.write('%.2f %i %.16f\n' % (iamp,avg,Nspikes_rms[i]))
    if varyE_bool==True:
        plt.errorbar(Is,Nspikes_avg,yerr=Nspikes_rms,capsize=2,label=r'E=%s' % str(elem)) # Flexible in case we want more decimals
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