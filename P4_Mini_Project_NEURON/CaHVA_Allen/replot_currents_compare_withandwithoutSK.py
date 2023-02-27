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
    
    gnaf   = 1.0
    gkaf   = 1.0
    gcahva = 0.2 # Available: [0.1,0.5,1.0,1.5,2.0]
    gsk    = 1.0
    gpas   = 1.0
    
    cm     = 1.0
    iamp   = 0.002 # Available: [0.002,0.005,0.01,0.015,0.02,0.04]  
    
    linewidths   = [1.0,1.5,1.5,1.5]
    mylinestyles = ['solid','dashed','dotted','dashdot']
    infolder1 = 'Shift_gCaHVA_V/Vshift_%i/Results/Soma%i/' % (vshift,somasize)+ 'current_idur%i_iamp'%idur+str(iamp)+'/'
    infolder2 = 'SK_Allen/Shift_gCaHVA_V/Vshift_%i/Results/Soma%i/' % (vshift,somasize)+ 'current_idur%i_iamp'%idur+str(iamp)+'/'
    folders = [infolder1,infolder2]
    legendstring = ['; No SK', '; With SK']
    Nfolders = len(folders)
      
    
    colors_1 = []
    colors_2 = []
    colors_3 = []
    colors_4 = []
    colors_5 = []
    colors_6 = []
    # 
    colors_1.append('tab:blue')
    colors_1.append('lightskyblue')
    colors_1.append('midnightblue')
    #
    colors_2.append('tab:orange')
    colors_2.append('moccasin')
    colors_2.append('darkgoldenrod')
    #
    colors_3.append('tab:red')
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
    
    plt.figure(figsize=(6,5))
    i = 0
    for i in range(Nfolders):
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
        if i>0:
            namestring = namestring + '_gSK'+str(gsk)+'p'    
        if vary_Ca_HVA==True:
            namestring = namestring + '_gCaHVA'+str(gcahva)+'p'    
        if vary_gpas==True: 
            namestring = namestring + '_gpas'+str(gpas)+'p'
        namestring = namestring +'_' 
        
        # Set names
        infolder = folders[i]
        filename = infolder+namestring+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_currents.txt'
        
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
        if i>0:
            I_sk   = data[:, 7]
        
    
        thiswidth     = linewidths[i]
        thislinestyle = mylinestyles[i]
        if plotnaf==True:
            plt.plot(time,I_naf,color=colors_1[i],label=r'$I_{\mathregular{naf}}$'+legendstring[i],linewidth=thiswidth,linestyle=thislinestyle)
        if plotkaf==True:
            plt.plot(time,I_kaf,color=colors_2[i],label=r'$I_{\mathregular{kaf}}$'+legendstring[i],linewidth=thiswidth,linestyle=thislinestyle)
        if plotcahva==True:
            plt.plot(time,I_CaHVA,color=colors_3[i],label=r'$I_{\mathregular{CaHVA}}$'+legendstring[i],linewidth=thiswidth,linestyle=thislinestyle)
        if plotcap==True:
            plt.plot(time,I_cap,color=colors_4[i],label=r'$I_{\mathregular{cap}}$'+legendstring[i],linewidth=thiswidth,linestyle=thislinestyle)
        if plotmem==True:
            plt.plot(time,I_mem,color=colors_5[i],label=r'$I_{\mathregular{mem}}$'+legendstring[i],linewidth=thiswidth,linestyle=thislinestyle)
        if plotpas==True:
            plt.plot(time,I_pas,color=colors_6[i],label=r'$I_{\mathregular{pas}}$'+legendstring[i],linewidth=thiswidth,linestyle=thislinestyle)
    plt.xlabel('Time [ms]')
    plt.ylabel('Current [nA]')
    plt.title(r'Currents, %s$\bar{g}_{\mathregular{CaHVA}}$, %s$\bar{g}_{\mathregular{SK}}$, $I_{\mathregular{input}}=$%s' % (str(gcahva),str(gsk),str(iamp)))
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.show()
