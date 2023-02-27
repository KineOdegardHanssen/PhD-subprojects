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

def manual(filename,idelay,idur,spikedurat):
    """Manual approach"""

    # Use numpy to read the trace data from the txt file
    data = numpy.loadtxt(filename)

    # Time is the first column
    time = data[:, 0]
    # Voltage is the second column
    voltage = data[:, 1]
    
    Vrest, Vrest_rms = avg_and_rms(voltage)
    
    return Vrest, Vrest_rms

if __name__ == '__main__':
    spikedurat = -40
    idur       = 1000 #100 # ms
    idelay     = 10
    v_init     = -86.8 # mV
    Ra         = 100
    somasize   = 10 # 15 # 
    dtexp      = -7
    
    cm = 1.0
    iamp = 0
    
    # Set names
    outfolder = 'Results/Soma%i/' % somasize
    outfilename_Vrest = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Vrest.txt'
    outfile_Vrest = open(outfilename_Vrest,'w')
    
    # Read file and get results
    infolder = outfolder + 'current_idur%i_iamp'%idur+str(iamp)+'/'
    filename = infolder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_V.txt' % dtexp
    Vrest, Vrest_rms = manual(filename,idelay,idur,spikedurat)
    
    # Write to file (Vrest_rms is a sort of safeguard)
    outfile_Vrest.write('%.16e %.16e' % (Vrest, Vrest_rms))
    
    outfile_Vrest.close()