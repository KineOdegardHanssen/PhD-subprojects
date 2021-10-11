import numpy as np
import matplotlib.pyplot as plt

factor = False
dtexp = -8
smallCm = True
cm = 1
model  = 478513407 # 488462965 # 485694403 # 478513437 #  [478513407,478513437,488462965] # 
varymech = 'NaV' #'SK'#'Im_v2'#'Kv2like'#'Kd'#'pas' #
varyE_bool = True
varyE = 30 # Choose one: [30,40,50,53,60,65,70]
varyg = 'None'
idur = 1000
iamp = 0.44
    
varylist = [] # Should be redundant
plotstring = '_vary'
if varyE_bool==True:
    varylist = varyE
    plotstring = plotstring + 'E'
else:
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
#testmodel = 488462965 # 478513437#478513407# 485694403 # 489931686 #  ### Weirdo:480633479#  
idelay   = 100.0
afteri   = 100.0
tstart   = -500.
tstop    = idur+afteri+idelay
v_init   = -86.5
somasize = 10 # 15 # 
dendlen  = 1000
denddiam = 2

currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
changestring  = ''
if varyE_bool==True:
    changestring = changestring+'_E'+str(varyE)+'_gdflt'
else:
    changestring = changestring+'_Edefault_g'+str(varyg)
infolder = 'Results/%i/IStim/Soma%i/dendlen%i/denddiam'% (model,somasize,dendlen)+str(denddiam)+'/'+ folderstring
infolder = infolder+currentfolder
if factor==True:
    infilename = infolder+'baspv_cmf'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_vinit'+str(v_init)+changestring+'_sprx_V.txt'
else:
    infilename = infolder+'baspv_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_vinit'+str(v_init)+changestring+'_sprx_V.txt'

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

plt.figure(figsize=(6,5))
if varyE_bool==True:
    plt.plot(t,V, label='E=%s, I=%.2f' % (str(varyE),iamp))
else:
    plt.plot(t,V, label='g=%s, I=%.2f' % (str(varyg),iamp))
plt.xlabel('time (ms)')
plt.ylabel('Soma potential (mV)')      
if varymech=='NaV':
    if varyE_bool==True:
        plt.title(r'fI-curve, perisomatic ball-and-stick, changing $E_{Na}$ soma and prox.dend.')
    else:
        plt.title(r'fI-curve, perisomatic ball-and-stick, changing $\bar{g}_{Na}$ soma and prox.dend.')
plt.legend(loc='upper left')
plt.show()
