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
    
    vmax = max(voltage) 
    vmin = min(voltage) 
    deltav = vmax-vmin
    vthr  = -20   # If there is a peak above this value, we count it # Old: vmax-0.15*deltav
    vprev = vthr-40 # A peak never kicks in at initiation, even if I change vthr
    durthr = spikedurat # Height at which we measure the duration. Need to be calibrated by peaks? # Was -20
    Npeaks = 0
    peakvals  = []
    peaktimes = []
    passtimes_up = [] # Hehe
    passvals_up  = []
    passtimes_down = []
    passvals_down  = []
    for i in range (1,len(voltage)-1):  
        if voltage[i-1]<voltage[i] and voltage[i+1]<voltage[i] and voltage[i]>vthr:
            peaktimes.append(time[i])
            peakvals.append(voltage[i])
            Npeaks+=1
        if voltage[i]>=durthr and voltage[i-1]<durthr: # Passing upwards
            tbef = time[i-1]
            taft = time[i]
            Vbef = voltage[i-1]
            Vaft = voltage[i]
            a = (Vaft-Vbef)/(taft-tbef)
            b = Vbef-a*tbef
            tint = (durthr-b)/a
            Vint = a*tint+b
            passtimes_up.append(tint)
            passvals_up.append(Vint) # For plotting
        elif voltage[i]>=durthr and voltage[i+1]<durthr: # Passing downwards
            tbef = time[i]
            taft = time[i+1]
            Vbef = voltage[i]
            Vaft = voltage[i+1]
            a = (Vaft-Vbef)/(taft-tbef)
            b = Vbef-a*tbef
            tint = (durthr-b)/a
            Vint = a*tint+b
            passtimes_down.append(tint)
            passvals_down.append(Vint) # For plotting
        if time[i]>idur+idelay:
            break
    
    # Checking if we've got consistent firing:
    if Npeaks!=0:
        if peaktimes[-1]<=(idur/2.+idelay): #Checking if there's no firing in the last half of the stim. interval
            Npeaks=0                        # Is that a proper limit? Last third instead?
    
    # I should probably plot stuff, to make sure... If-test?
    dur = []
    isi = []
    Ndur = min([len(passtimes_up),len(passtimes_down)])
    for i in range(Ndur):
        dur.append(passtimes_down[i]-passtimes_up[i])
    for i in range(1,len(peaktimes)):
        isi.append(peaktimes[i]-peaktimes[i-1])
    isi_first = isi[:5]
    isi_last = isi[-5:]
    
    ## Avg and rms:
    isi_first_avg, isi_first_rms = avg_and_rms(isi_first)
    isi_last_avg, isi_last_rms = avg_and_rms(isi_last)
    
    return Npeaks, isi_first_avg, isi_first_rms, isi_last_avg, isi_last_rms

if __name__ == '__main__':
    spikedurat = -40
    idur       = 1000 #100 # ms
    idelay     = 10
    v_init     = -86 # mV
    Ra         = 100
    somasize   = 10 # 15 # 
    dtexp      = -7
    t_before_rec = -600.
    
    model_folder = '' # We are in the model folder
    
    cm = 1.0
    iamps = [0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03,0.031,0.032,0.033,0.034,0.035,0.036,0.037,0.038,0.039,0.04,0.041,0.042,0.043,0.044,0.045,0.046,0.047,0.048,0.049,0.05,0.051,0.052,0.053,0.054,0.055,0.056,0.057,0.058,0.059,0.06]#[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5]
    gcahvas  = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2.0,2.5]
    
    for gcahva in gcahvas:
        namestring = '_gCaHVA'+str(gcahva)+'p_'
        
        Namps = len(iamps)
        
        Nspikes = numpy.zeros(Namps)
        isi_ratio_avg = numpy.zeros(Namps)
        isi_ratio_rms = numpy.zeros(Namps)
        isi_first_avg = numpy.zeros(Namps)
        isi_first_rms = numpy.zeros(Namps)
        isi_last_avg = numpy.zeros(Namps)
        isi_last_rms = numpy.zeros(Namps)
        
        # Set names
        outfolder = 'Results/Soma%i/' % somasize
        outfilename = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+namestring+'_manual_cm'+str(cm)+'_ISIratio_vs_I.txt'
        plotname_ratio   = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+namestring+'_manual_cm'+str(cm)+'_ISIratio_vs_I.png'
        plotname_firstlast = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+namestring+'_manual_cm'+str(cm)+'_ISIfirstlast_vs_I.png'
        # make files
        outfile = open(outfilename,'w')
        firedbefore = False
        for j in range(Namps):
            iamp = iamps[j]
            print('Step ', j+1, ' of', Namps, 'cm:', cm, 'iamp:', iamp)
            infolder = outfolder + 'current_idur%i_iamp'%idur+str(iamp)+'/'
            filename = infolder+namestring+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_V.txt'
            #try:
            #print('In try')
            Nspikes[j], isi_first_avg[j], isi_first_rms[j], isi_last_avg[j], isi_last_rms[j] = manual(filename,idelay,idur,spikedurat)
            if Nspikes[j]!=0 and firedbefore==False:
                firedbefore = True
            if Nspikes[j]!=0:
                isi_ratio_avg[j] = isi_last_avg[j]/isi_first_avg[j]
                isi_ratio_rms[j] = numpy.sqrt(isi_ratio_avg[j]**2*((isi_first_rms[j]/isi_first_avg[j])**2+(isi_last_rms[j]/isi_last_avg[j])**2))
                outfile.write('%.5f %.10f %.10f %.10f %.10f %.10f %.10f\n' % (iamp,isi_ratio_avg[j],isi_ratio_rms[j],isi_first_avg[j],isi_first_rms[j],isi_last_avg[j],isi_last_rms[j]))
        outfile.close()
        
        # Plot results
        plt.figure(figsize=(6,5))
        plt.errorbar(iamps,isi_ratio_avg, yerr=isi_ratio_rms, capsize=2)
        plt.xlabel(r'$I$ [nA]')
        plt.ylabel(r'Ratio')
        plt.title(r'Ratio $\frac{\mathregular{ISI}_{\mathregular{last}}}{\mathregular{ISI}_{\mathregular{first}}}$')
        plt.tight_layout()    
        plt.savefig(plotname_ratio)
        
        plt.figure(figsize=(6,5))
        plt.errorbar(iamps,isi_first_avg, yerr=isi_first_rms, capsize=2,label='first')
        plt.errorbar(iamps,isi_last_avg, yerr=isi_last_rms, capsize=2,label='last')
        plt.xlabel(r'$I$ [nA]')    
        plt.ylabel(r'ISI [ms]')
        plt.title(r'$\mathregular{ISI}_{\mathregular{first}}$ and $\mathregular{ISI}_{\mathregular{last}}$')
        plt.legend(loc='upper right')
        plt.tight_layout()
        plt.savefig(plotname_firstlast)
