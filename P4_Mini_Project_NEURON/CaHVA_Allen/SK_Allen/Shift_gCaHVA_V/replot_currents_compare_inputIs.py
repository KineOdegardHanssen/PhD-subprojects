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
    I_SK    = data[:, 7]
    #t[k],I_naf[k],I_kaf[k],I_Ca_HVA[k],icap[k],imem[k],ipas[k]
    return time, I_naf, I_kaf, I_CaHVA, I_cap, I_mem, I_pas, I_SK
    
    

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
    vary_SK      = True # False
    vary_Ca_HVA  = True # False
    vary_gpas    = False # 

    ### Plotting options:
    plotnaf   = False # True #
    plotkaf   = False # True #
    plotsk    = False # True  # 
    plotcahva = True # False # 
    plotcap   = False # True #
    plotmem   = False # True #
    plotpas   = False # True #
    
    gnaf    = 1.0
    gkaf    = 1.0
    gcahva  = 0.2
    gsk     = 0.5 # Available: [0.1,0.5,1.0,1.5,2.0]
    gpas    = 1.0
    
    linewidths   = [1.0,1.5,1.5,1.5]
    mylinestyles = ['solid','solid','solid','solid']#['solid','dashed','dotted','dashdot']
    
    cm    = 1.0
    iamps = [0.002,0.01] # Available: [0.002,0.005,0.01,0.015,0.02,0.04]    
    
    I_naf   = []
    I_kaf   = []
    I_CaHVA = []
    I_cap   = []
    I_mem   = []
    I_pas   = []
    I_SK    = []
    
    colors_1 = []
    colors_2 = []
    colors_3 = []
    colors_4 = []
    colors_5 = []
    colors_6 = []
    colors_7 = []
    # First
    colors_1.append('tab:blue')
    colors_1.append('lightskyblue')
    colors_1.append('midnightblue')
    #
    colors_2.append('tab:orange')
    colors_2.append('moccasin')
    colors_2.append('darkgoldenrod')
    #
    colors_3.append('tab:red')
    colors_3.append('lightcoral')
    colors_3.append('maroon')
    #
    colors_4.append('tab:green')
    colors_4.append('lightgreen')
    colors_4.append('darkgreen')
    #
    colors_5.append('tab:brown')
    colors_5.append('peru')
    colors_5.append('saddlebrown')
    #
    colors_6.append('tab:purple')
    colors_6.append('orchid')
    colors_6.append('purple')
    #
    colors_7.append('tab:gray')
    colors_7.append('lightgray')
    colors_7.append('black')
    
    plt.figure(figsize=(6,5))
    i = 0
    for iamp in iamps:
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
        if vary_SK==True:
            namestring = namestring + '_gSK'+str(gsk)+'p'  
        if vary_Ca_HVA==True:
            namestring = namestring + '_gCaHVA'+str(gcahva)+'p'    
        if vary_gpas==True: 
            namestring = namestring + '_gpas'+str(gpas)+'p'
        namestring = namestring +'_' 
        
        # Set names
        infolder = 'Vshift_%i/Results/Soma%i/' % (vshift,somasize)+ 'current_idur%i_iamp'%idur+str(iamp)+'/'
        
        filename = infolder+namestring+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_currents.txt'
        time,I_naf_this,I_kaf_this,I_CaHVA_this,I_cap_this,I_mem_this,I_pas_this, I_SK_this = manual(filename)
        
        I_naf.append(I_naf_this)
        I_kaf.append(I_kaf_this)
        I_CaHVA.append(I_CaHVA_this)
        I_cap.append(I_cap_this)
        I_mem.append(I_mem_this)
        I_pas.append(I_pas_this)
        I_SK.append(I_SK_this)
        
        thiswidth     = linewidths[i]
        thislinestyle = mylinestyles[i]
        if plotnaf==True:
            plt.plot(time,I_naf_this,color=colors_1[i],label=r'$I_{\mathregular{naf}}$, $I_{\mathregular{input}}=$%s' % str(iamp),linewidth=thiswidth,linestyle=thislinestyle)
        if plotkaf==True:
            plt.plot(time,I_kaf_this,color=colors_2[i],label=r'$I_{\mathregular{kaf}}$, $I_{\mathregular{input}}=$%s' % str(iamp),linewidth=thiswidth,linestyle=thislinestyle)
        if plotsk==True:
            plt.plot(time,I_SK_this,color=colors_7[i],label=r'$I_{\mathregular{SK}}$, $I_{\mathregular{input}}=$%s' % str(iamp),linewidth=thiswidth,linestyle=thislinestyle)
        if plotcahva==True:
            plt.plot(time,I_CaHVA_this,color=colors_3[i],label=r'$I_{\mathregular{CaHVA}}$, $I_{\mathregular{input}}=$%s' % str(iamp),linewidth=thiswidth,linestyle=thislinestyle)
        if plotcap==True:
            plt.plot(time,I_cap_this,color=colors_4[i],label=r'$I_{\mathregular{cap}}$, $I_{\mathregular{input}}=$%s' % str(iamp),linewidth=thiswidth,linestyle=thislinestyle)
        if plotmem==True:
            plt.plot(time,I_mem_this,color=colors_5[i],label=r'$I_{\mathregular{mem}}$, $I_{\mathregular{input}}=$%s' % str(iamp),linewidth=thiswidth,linestyle=thislinestyle)
        if plotpas==True:
            plt.plot(time,I_pas_this,color=colors_6[i],label=r'$I_{\mathregular{pas}}$, $I_{\mathregular{input}}=$%s' % str(iamp),linewidth=thiswidth,linestyle=thislinestyle)
        i+=1
    plt.xlabel('Time [ms]')
    plt.ylabel('Current [nA]')
    plt.title('Currents')
    plt.legend(loc='lower right')#'upper left')
    plt.tight_layout()
    plt.show()
