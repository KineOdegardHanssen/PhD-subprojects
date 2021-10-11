import numpy as np
import matplotlib.pyplot as plt

############### Somaonly-HH ###########################################################
cms        = [1.0]
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

Nspikes = []
for i in range(Ncm):
    cm = cms[i]
    theseNspikes = []
    for j in range(Ni):
        iamp = iamps[j]
        folder = 'Somaonly/Results/IStim/Soma%i/current_idur%i_iamp'%(somasize,idur)+str(iamp)+'/'
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
    plt.plot(iamps,theseNspikes,'-o',label=r'%s*$C_m$' % (str(cm)))
    Nspikes.append(theseNspikes)

fIsomaonlyHH = theseNspikes

############### Somaonly-PV ###########################################################
testmodels = [488462965]#[478513437,478513407,478513437,488462965]
cm_factors = [1]
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
        folder = 'Somaonly/Results/IStim/Soma%i/current_idur%i_iamp'%(somasize,idur)+str(iamp)+'/'
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

fIsomaonlyPV_avg = Nspikes_avg
fIsomaonlyPV_rms = Nspikes_rms

############### BAS-HH ################################################################
dtexp = -8
smallCm = True
cm = 1
varymech = 'Na' #'K'#'leak'
varyE_bool = True
varyE = [50] # Default
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

Is = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44] # Need to hard-code :( # Should I have more Is, though?
NI = len(Is)
cm        = 1.0 
cmfac     = 1.0

outfolder = 'Ball-and-stick models/BAS_somaHH_dendHH/Results/IStim/Soma%i/' % somasize

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
        infolder = 'Ball-and-stick models/BAS_somaHH_dendHH/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+ folderstring
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

fIBASHH = Nspikes

############### BAS-PV ################################################################
smallCm = True
cm = 1
varymech = 'Kd'#'NaV' #'SK'#'Im_v2'#'Kv2like'#'pas' #
varyE_bool = True
varyE = [-107]
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
models  = [489931686] # [478513407,478513437,489931686,485694403,488462965]
Nmodels = len(models)

outfolder = 'Ball-and-stick models/BAS_somaPV_dendPV/Results/Soma%i/' % somasize
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
            infolder = 'Ball-and-stick models/BAS_somaPV_dendPV/Results/%i/IStim/Soma%i/dendlen%i/denddiam'% (model,somasize,dendlen)+str(denddiam)+'/'+ folderstring
            infolder = infolder+currentfolder
            infilename = infolder+'baspv_cmf'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_vinit'+str(v_init)+changestring+'_sprx_V.txt'
            outfilename = infolder+'baspv_cmf'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_vinit'+str(v_init)+changestring+'_sprx_If.txt'
            
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
        Nspikes_rms[i] = np.sqrt(rms/(Nmodels-1))
        outfile_all.write('%.2f %i %.16f\n' % (iamp,avg,Nspikes_rms[i]))

fIBASPV_avg = Nspikes_avg
fIBASPV_rms = Nspikes_rms

############### Allen ################################################################
dtexp = -8
smallCm = True
cmfs = [1.0]
idur = 1000
idelay   = 100.0
afteri   = 100.0
tstart   = -500.
tstop    = idur+afteri+idelay
v_init   = -86.5
somasize = 10
dendlen  = 1000
denddiam = 2

NI   = 23 # 10 # 
Imax = 0.18
if NI==23:
    Imax = 0.44
Is   = np.linspace(0,Imax,NI)
models  = [478513407,478513437,488462965] # [478513437,488462965] # 
Nmodels = len(models)
Ncms    = len(cms)

plt.figure(figsize=(6,5))
for l in range(Ncms):
    cmf = cmfs[l]
    Nspikes_avg = np.zeros(NI)
    Nspikes_rms = np.zeros(NI)
    for i in range(NI):
        iamp = Is[i]
        if iamp==0.0:
            iamp=0
        #if cmf==1:
        #    cmf = 1.0
        theseNspikes = []
        currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
        for j in range(Nmodels):
            testmodel = models[j]
            if testmodel==480633479:
                v_init = -96.8#-83.7#-90#-86.5# # Have a test here too
                cm_soma = 0.704866 # 0.704866118957
                cm_dend = 0.704866 # 0.704866118957
                cm_axon = 0.704866 # 0.704866118957
            elif testmodel==496497595:
                v_init = -86.5
                cm_soma = 1.14805
                cm_dend = 9.98231
                cm_axon = 3.00603
            elif testmodel==488462965:
                v_init = -86.5 # Maybe I should have changed this...
                cm_soma = 3.31732779736 # Strange values, right?
                cm_dend = 3.31732779736
                cm_axon = 3.31732779736
            elif testmodel==497230463:
                v_init = -90
                cm_soma = 1.23729
                cm_dend = 2.57923
                cm_axon = 5.02697
            elif testmodel==497233075:
                v_init = -90
                cm_soma = 1.64168
                cm_dend = 2.83035
                cm_axon = 9.98442
            elif testmodel==478513437:
                v_init = -86.8
                cm_soma = 2.34539964752
                cm_dend = 2.34539964752
                cm_axon = 2.34539964752
            elif testmodel==478513407:
                v_init = -83.7
                cm_soma = 1.0
                cm_dend = 1.0
                cm_axon = 1.0
            elif testmodel==497233271:
                v_init = -90
                cm_soma = 0.783229
                cm_dend = 1.94512
                cm_axon = 8.25387
            elif testmodel==489931686:
                v_init = -95.7
                cm_soma = 1.66244903951
                cm_dend = 1.66244903951
                cm_axon = 1.66244903951
            infilename = "Allen_test_changecapacitance/figures/%i/current_idur%i_iamp" % (testmodel,idur) + str(iamp)+"/idur%i_iamp" % idur + str(iamp)+"_changecmf" + str(cmf) + "_cmdend" + str(cm_dend) + "_cmaxon"+ str(cm_axon) + "_vinit"+str(v_init)+"_addedRa.txt"
            
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
            print(Npeaks, ' peaks for model ', testmodel, ', current ', iamp)
            print('vmax:', vmax)
            
            Nspikes_avg[i]+=Npeaks
            theseNspikes.append(Npeaks)
        Nspikes_avg[i]/=Nmodels
        avg = Nspikes_avg[i]
        rms = 0
        for k in range(Nmodels):
            rms+=(avg-theseNspikes[k])**2
        Nspikes_rms[i] = np.sqrt(rms/(Nmodels-1))

fIAllen_avg = Nspikes_avg
fIAllen_rms = Nspikes_rms


############### Plotting ################################################################

#fIsomaonlyHH
##
#fIsomaonlyPV_avg
#fIsomaonlyPV_rms
##
#fBASHH = Nspikes
##
#fIBASPV_avg 
#fIBASPV_rms 
##
#fIAllen_avg 
#fIAllen_rms

plotname = 'Comparemodels/all_cmf1_standard_l1000_d2_fI.png'

print('len(iamps):',len(iamps))
print('len(fIAllen_avg):',len(fIAllen_avg))
print('len(fIAllen_rms):', len(fIAllen_rms))

plt.figure(figsize=(6,5))
plt.plot(iamps,fIsomaonlyHH,label='One comp., HH')
#plt.errorbar(iamps,fIsomaonlyPV_avg, yerr=fIsomaonlyPV_rms,capsize=2,label='One comp., PV')
plt.errorbar(iamps,fIsomaonlyPV_avg, yerr=fIsomaonlyPV_rms,capsize=2,label='One comp., PV')
plt.plot(iamps,fIBASHH,label='Ball-and-stick, HH')
plt.errorbar(iamps,fIBASPV_avg, yerr=fIBASPV_rms,capsize=2,label='Ball-and-stick, PV')
plt.errorbar(iamps,fIAllen_avg, yerr=fIAllen_rms,capsize=2,label='Biophys.')
plt.xlabel('I (nA)')
plt.ylabel('f (Hz)')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig(plotname)
