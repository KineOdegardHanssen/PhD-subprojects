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
    plt.legend(loc='lower left')
    plt.tight_layout()
    plt.show() 

if __name__ == '__main__':
    spikedurat = -40
    idur       = 1000 #100 # ms
    idelay     = 10
    v_init     = -86.8 # mV
    Ra         = 100
    somasize   = 10 # 15 # 
    dtexp      = -7
    
    outfolder = 'Results/Soma%i/' % somasize # We are in the model folder
    
    cm = 1.0
    iamps = [0.051]#[0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03,0.031,0.032,0.033,0.034,0.035,0.036,0.037,0.038,0.039,0.04,0.041,0.042,0.043,0.044,0.045,0.046,0.047,0.048,0.049,0.05,0.051,0.052,0.053,0.054,0.055,0.056,0.057,0.058,0.059,0.06]#
    
    Namps = len(iamps)
    
    for j in range(Namps):
        iamp = iamps[j]
        print('Step ', j+1, ' of', Namps, 'cm:', cm, 'iamp:', iamp)
        infolder = outfolder + 'current_idur%i_iamp'%idur+str(iamp)+'/'
        filename = infolder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_V.txt' % dtexp
        replot(filename,iamp,idelay,idur,spikedurat)
