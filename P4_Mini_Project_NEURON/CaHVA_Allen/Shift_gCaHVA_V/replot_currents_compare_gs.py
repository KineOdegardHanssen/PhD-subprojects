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

def manual(filename):
    """Manual approach"""

    # Use numpy to read the trace data from the txt file
    data = numpy.loadtxt(filename)

    # Time is the first column
    time = data[:, 0]
    # Currents in the other columns
    I_naf   = data[:, 1]
    I_kaf   = data[:, 2]
    I_CaHVA = data[:, 3]
    I_cap   = data[:, 4]
    I_mem   = data[:, 5]
    I_pas   = data[:, 6]
    #t[k],I_naf[k],I_kaf[k],I_Ca_HVA[k],icap[k],imem[k],ipas[k]
    return time, I_naf, I_kaf, I_CaHVA, I_cap, I_mem, I_pas
    
    

if __name__ == '__main__':
    vshift     = 0
    skiptime   = 500 # 200
    spikedurat = -40
    idur       = 1000 #100 # ms
    idelay     = 100
    v_init     = -86 # mV
    Ra         = 100
    somasize   = 10 # 15 # 
    dtexp      = -7
    t_before_rec = -600.
    
    
    model_folder = '' # We are in the model folder
    
    vary_naf     = False
    vary_kaf     = False
    vary_Ca_HVA  = True # False
    vary_gpas    = False # 

    ### Plotting options:
    plotnaf   = False # True # 
    plotkaf   = False # True # 
    plotcahva = True
    plotcap   = False # True # 
    plotmem   = False # True # 
    plotpas   = False # True # 
    
    gnaf    = 1.0
    gkaf    = 1.0
    gcahvas = [0.1,0.5] # Available: [0.1,0.5,1.0,1.5,2.0]
    gpas    = 1.0
    
    linewidths   = [1.0,1.5,1.5,1.5]
    mylinestyles = ['solid','dashed','dotted','dashdot']
    
    cm    = 1.0
    iamp  = 0.002 # Available: [0.002,0.005,0.01,0.015,0.02,0.04]    
    
    I_naf   = []
    I_kaf   = []
    I_CaHVA = []
    I_cap   = []
    I_mem   = []
    I_pas   = []
    
    plt.figure(figsize=(6,5))
    i = 0
    for gcahva in gcahvas:
        varyE = 0
        varymech = 'None' # 'Epas' # 'EK' # 'ENa' # 
        namestring = ''
        if varymech=='ENa':
            varyE = 63 #[40,50,60,70]
            namestring = namestring + 'ENa'+str(varyE)
        elif varymech=='EK':
            varyE = -97#-107
            namestring = namestring + 'EK'+str(varyE)
        elif varymech=='Epas':
            varyE = -20 # Vary by shifts
            namestring = namestring + 'Epasshift'+str(varyE)
        if vary_naf==True:
            namestring = namestring + '_gnaf'+str(gnaf)+'p'
        if vary_kaf==True:
            namestring = namestring + '_gkaf'+str(gkaf)+'p'
        if vary_Ca_HVA==True:
            namestring = namestring + '_gCaHVA'+str(gcahva)+'p'    
        if vary_gpas==True: 
            namestring = namestring + '_gpas'+str(gpas)+'p'
        namestring = namestring +'_' 
        
        # Set names
        infolder = 'Vshift_%i/Results/Soma%i/' % (vshift,somasize)+ 'current_idur%i_iamp'%idur+str(iamp)+'/'
        
        filename = infolder+namestring+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_currents.txt'
        time,I_naf_this,I_kaf_this,I_CaHVA_this,I_cap_this,I_mem_this,I_pas_this = manual(filename)
        
        I_naf.append(I_naf_this)
        I_kaf.append(I_kaf_this)
        I_CaHVA.append(I_CaHVA_this)
        I_cap.append(I_cap_this)
        I_mem.append(I_mem_this)
        I_pas.append(I_pas_this)
    
        thiswidth     = linewidths[i]
        thislinestyle = mylinestyles[i]
        if plotnaf==True:
            plt.plot(time,I_naf_this,color='tab:blue',label=r'$I_{\mathregular{naf}}$, %s$\bar{g}_{\mathregular{CaHVA}}$' % str(gcahva),linewidth=thiswidth,linestyle=thislinestyle)
        if plotkaf==True:
            plt.plot(time,I_kaf_this,color='tab:orange',label=r'$I_{\mathregular{kaf}}$, %s$\bar{g}_{\mathregular{CaHVA}}$' % str(gcahva),linewidth=thiswidth,linestyle=thislinestyle)
        if plotcahva==True:
            plt.plot(time,I_CaHVA_this,color='tab:red',label=r'$I_{\mathregular{CaHVA}}$, %s$\bar{g}_{\mathregular{CaHVA}}$' % str(gcahva),linewidth=thiswidth,linestyle=thislinestyle)
        if plotcap==True:
            plt.plot(time,I_cap_this,color='tab:green',label=r'$I_{\mathregular{cap}}$, %s$\bar{g}_{\mathregular{CaHVA}}$' % str(gcahva),linewidth=thiswidth,linestyle=thislinestyle)
        if plotmem==True:
            plt.plot(time,I_mem_this,color='tab:olive',label=r'$I_{\mathregular{mem}}$, %s$\bar{g}_{\mathregular{CaHVA}}$' % str(gcahva),linewidth=thiswidth,linestyle=thislinestyle)
        if plotpas==True:
            plt.plot(time,I_pas_this,color='tab:purple',label=r'$I_{\mathregular{pas}}$, %s$\bar{g}_{\mathregular{CaHVA}}$' % str(gcahva),linewidth=thiswidth,linestyle=thislinestyle)
        i+=1
    plt.xlabel('Time [ms]')
    plt.ylabel('Current [nA]')
    plt.title('Currents')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.show()
