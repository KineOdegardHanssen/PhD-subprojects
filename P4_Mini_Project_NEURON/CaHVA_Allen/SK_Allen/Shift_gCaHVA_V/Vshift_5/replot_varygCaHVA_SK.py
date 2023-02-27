import numpy
import matplotlib.pyplot as plt

def replot(filename,iamp,idelay,idur,spikedurat):
    """Manual approach"""

    # Use numpy to read the trace data from the txt file
    data = numpy.loadtxt(filename)

    # Time is the first column
    time = data[:, 0]
    # Voltage is the second column
    voltage = data[:, 1]
    
    plt.figure(figsize=(6,5))
    plt.plot(time,voltage)
    plt.xlabel('Time [ms]')
    plt.ylabel('Voltage [mV]')
    plt.title('I=%s nA' % str(iamp))
    plt.tight_layout()
    plt.show() 

if __name__ == '__main__':
    spikedurat = -40
    idur       = 1000 #100 # ms
    idelay     = 10
    v_init     = -86 # mV
    Ra         = 100
    somasize   = 10 # 15 # 
    dtexp      = -7
    t_before_rec = -600.
    
    vary_naf     = False
    vary_kaf     = False
    vary_SK      = True # False # 
    vary_Ca_HVA  = True # False
    vary_gpas    = False # 
    
    gnaf   = 1.0
    gkaf   = 1.0
    gcahva = 0.1
    gsk    = 1.0
    gpas   = 1.0
    
    outfolder = 'Results/Soma%i/' % somasize
    
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
    namestring = namestring + '_'
    
    cm = 1.0
    iamps = [0.051]#[0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03,0.031,0.032,0.033,0.034,0.035,0.036,0.037,0.038,0.039,0.04,0.041,0.042,0.043,0.044,0.045,0.046,0.047,0.048,0.049,0.05,0.051,0.052,0.053,0.054,0.055,0.056,0.057,0.058,0.059,0.06]#
    
    Namps = len(iamps)
    
    for j in range(Namps):
        iamp = iamps[j]
        print('Step ', j+1, ' of', Namps, 'cm:', cm, 'iamp:', iamp)
        infolder = outfolder + 'current_idur%i_iamp'%idur+str(iamp)+'/'
        filename = infolder+namestring+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_V.txt'
        replot(filename,iamp,idelay,idur,spikedurat)
