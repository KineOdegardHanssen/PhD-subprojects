import numpy as np
import matplotlib.pyplot as plt

testmodel = 478513407#488462965#478513437#
varymech = 'NaV'
idur = 1000
somasize = 10
cm = 1.0

if varymech=='NaV':
    folderstring = 'VaryNa/' 
    potentials = [-60,-53,-30,-10,10,30,53,60]
    g = 0.9 # default 0.0409177
    potstring ='Ena'
    textsnippet = '_varyNa'
elif varymech=='pas':
    folderstring = 'VaryPas/'
    textsnippet = '_varypas'
elif varymech=='Kd':
    folderstring = 'VaryKd/'
    textsnippet = '_varyKd'
elif varymech=='Kv2like':
    folderstring = 'VaryKv2like/'
    textsnippet = '_varyKv2like'
elif varymech=='Kv3_1':
    folderstring = 'VaryKv3_1/'
    textsnippet = '_varyKv3_1'
elif varymech=='SK':
    folderstring = 'VarySK/'
    textsnippet = '_varySK'
    potentials = [-65]
elif varymech=='K_T':
    folderstring = 'VaryK_T/'
    textsnippet = '_varyK_T'
elif varymech=='Im_v2':
    folderstring = 'VaryIm_v2/'
    textsnippet = '_varyIm_v2'


folder_base = 'Results/%i/IStim/' % testmodel + folderstring


plotname = folder_base+'baspv_idur%i_cm'% idur+str(cm)+'_varyparam_%s_g' % potstring+str(g)+'_Nspikes_vs_iamp.png'

plt.figure(figsize=(6,5))
for potential in potentials:
    infilename = folder_base+'basps_cellmodel%i_current_idur%i'% (testmodel,idur)+textsnippet+'_E'+str(potential)+'_g'+str(g)+'_Nspikes_vs_iamp_compare.txt'
    
    infile = open(infilename,'r')
    lines = infile.readlines()
    
    iamp = []
    Nspikes = []
    
    for line in lines:
        words = line.split()
        if len(words)>0:
            iamp.append(float(words[0]))
            Nspikes.append(float(words[1]))
    infile.close()
    
    plt.plot(iamp, Nspikes, label=r'$E$=%s' % str(potential))
    #plt.plot(iamp, Nspikes, '-o', label=r'$E$=%s' % str(potential))

plt.xlabel(r'$I$ (nA)')
plt.ylabel('Number of spikes (frequency)')
if  varymech=='NaV':
    plt.title(r'Number of spikes vs $E_{Na}$')
elif varymech=='SK' or varymech=='Kd' or varymech=='Kv2like' or varymech=='Kv3_1' or varymech=='K_T' or varymech=='Im_v2':
    plt.title(r'Number of spikes vs $E_K$')
elif varymech=='pas':
    plt.title(r'Number of spikes vs $E_l$')
plt.legend(loc='upper left')
plt.savefig(plotname)
plt.show()