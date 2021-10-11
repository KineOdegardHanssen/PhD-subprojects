"""Modified from 'Basic example 1 for eFEL'."""

import efel
import numpy
import matplotlib.pyplot as plt

def avg_and_rms(x):
    N = len(x)
    avgx = numpy.mean(x)
    rmsx = 0
    for i in range(N):
        rmsx += (avgx-x[i])**2
    rmsx = numpy.sqrt(rmsx/(N-1))
    return avgx,rmsx

    
    
if __name__ == '__main__':
    testmodel = 489931686 # 485694403 # 488462965 # 478513407 # 478513437 # 
    cm   = [0.8,0.9]
    NCms = len(cm)
    
    idur      = 1000 #100 # ms
    idelay    = 100
    iamp      = 0.4 # nA
    v_init    = -65 # mV
    Ra        = 100
    somasize  = 10 # 15 # 
    dendlen   = 1000
    denddiam  = 2
    nsegments = 200 
    
    varymech = 'Kd'#'NaV' #'SK'#'Im_v2'#'Kv2like'#'pas' #
    varyE_bool = True
    varyE = -107 #[-90,-80]#[-130,-107,-100,-70]#[40,50,53,65,70]#[-60,-53,-30,-10,10,30,60] # Default...407Kd-107 Default478513407pas: -83.6528;  Default478513407Na: 53
    varyg = 'None' # Default...407Kd2.94396e-010 # Default478513407pas: 0.000362109; Default478513407Na: 0.0409177
    
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
        plotstring   = plotstring + '_NaV'
    elif varymech=='pas':
        folderstring = 'VaryPas/'
        plotstring   = plotstring + '_Pas'
    elif varymech=='Kd':
        folderstring = 'VaryKd/'
        plotstring   = plotstring + '_Kd'
    elif varymech=='Kv2like':
        folderstring = 'VaryKv2like/'
        plotstring   = plotstring + '_Kv2like'
    elif varymech=='Kv3_1':
        folderstring = 'VaryKv3_1/'
        plotstring   = plotstring + '_Kv3_1'
    elif varymech=='SK':
        folderstring = 'VarySK/'
        plotstring   = plotstring + '_SK'
    elif varymech=='K_T':
        folderstring = 'VaryK_T/'
        plotstring   = plotstring + '_K_T'
    elif varymech=='Im_v2':
        folderstring = 'VaryIm_v2/'
        plotstring   = plotstring + '_Im_v2'

    changestring =''
    if varyE_bool==True:
        changestring = changestring+'_E'+str(varyE)+'_gdf'
    else:
        changestring = changestring+'_Edf_g'+str(varyg)
    
    
    currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    for j in range(NCms):
        print('Step ', j+1, ' of', NCms)
        infolder = 'Results/%i/IStim/Soma%i/dendlen%i/denddiam'% (testmodel,somasize,dendlen)+str(denddiam)+'/'+ folderstring+currentfolder
        filename = infolder+'baspv_cmf'+str(cm[j])+'_idur%i_iamp'%idur+str(iamp)+changestring+'_sprx_V.txt' 
        
        data = numpy.loadtxt(filename)

        # Time is the first column
        time = data[:, 0]
        # Voltage is the second column
        voltage = data[:, 1]
        
        plt.figure(figsize=(6,5))        
        plt.plot(time,voltage)
        plt.xlabel('Time (ms)')
        plt.ylabel('Voltage (mV)')
        plt.title('Time vs voltage, $C_m$=%s of soma and proximal dendrite' % (str(cm[j])))
        plt.show()