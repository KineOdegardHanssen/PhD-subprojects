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

def manual(filename,iamp,idelay,idur,spikedurat,skiptime,gcahva,plotnaf,plotkaf,plotcahva,plotcap,plotmem,plotpas):
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
    
    plt.figure(figsize=(6,5))
    if plotnaf==True:
        plt.plot(time,I_naf,label=r'$I_{\mathregular{naf}}$')
    if plotkaf==True:
        plt.plot(time,I_kaf,label=r'$I_{\mathregular{kaf}}$')
    if plotcahva==True:
        plt.plot(time,I_CaHVA,label=r'$I_{\mathregular{CaHVA}}$')
    if plotcap==True:
        plt.plot(time,I_cap,label=r'$I_{\mathregular{cap}}$')
    if plotmem==True:
        plt.plot(time,I_mem,label=r'$I_{\mathregular{mem}}$')
    if plotpas==True:
        plt.plot(time,I_pas,label=r'$I_{\mathregular{pas}}$')
    plt.xlabel('Time [ms]')
    plt.ylabel('Current [nA]')
    plt.title('Currents')
    plt.legend(loc='lower left')
    plt.tight_layout()
    plt.show()
    

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
    plotnaf   = True
    plotkaf   = True
    plotcahva = True
    plotcap   = True
    plotmem   = True
    plotpas   = True
    
    gnaf    = 1.0
    gkaf    = 1.0
    gcahvas = [0.1,0.5]
    gpas    = 1.0
    
    cm = 1.0
    iamps =  [0.002] # Available: [0.002,0.005,0.01,0.015,0.02,0.04]    
    Namps = len(iamps)
    
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
        outfolder = 'Vshift_%i/Results/Soma%i/' % (vshift,somasize)
        for j in range(Namps):
            iamp = iamps[j]
            print('Step ', j+1, ' of', Namps, 'cm:', cm, 'iamp:', iamp)
            infolder = outfolder + 'current_idur%i_iamp'%idur+str(iamp)+'/'
            filename = infolder+namestring+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_currents.txt'
            manual(filename,iamp,idelay,idur,spikedurat,skiptime,gcahva,plotnaf,plotkaf,plotcahva,plotcap,plotmem,plotpas)
            
