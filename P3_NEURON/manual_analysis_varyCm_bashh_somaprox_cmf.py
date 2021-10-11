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
    
    vmax = max(voltage) 
    vmin = min(voltage) 
    deltav = vmax-vmin
    vthr  = -20   # If there is a peak above this value, we count it # Old: vmax-0.15*deltav
    vprev = vthr-40 # A peak never kicks in at initiation, even if I change vthr
    durthr = spikedurat # Height at which we measure the duration. Need to be calibrated by peaks? # Was -20
    Npeaks = 0
    peakvals = []
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
    
    # Checking if we've got consistent firing:
    if peaktimes[-1]<=(idur/2+idelay): # Checking if there's no firing in the last half of the stim. interval
        Npeaks=0                       # Is that a proper limit? Last third instead?
    
    # I should probably plot stuff, to make sure... If-test?
    dur = []
    for i in range(len(passtimes_up)):
        dur.append(passtimes_down[i]-passtimes_up[i])
    
    '''
    plt.figure(figsize=(6,5))
    plt.plot(time,voltage,',')
    plt.plot(peaktimes,peakvals,'o',label='peaks')
    plt.plot(passtimes_up,passvals_up,'o',label='dur basis, up')
    plt.plot(passtimes_down,passvals_down,'o',label='dur basis, down')
    plt.xlabel('Time [ms]')
    plt.ylabel('Voltage [mV]')
    plt.title('Testing implementation')
    plt.legend(loc='lower left')
    plt.tight_layout()
    plt.show() 
    '''
    
    print('dur:',dur)
    
    ## Avg and rms:
    peakvals_avg, peakvals_rms = avg_and_rms(peakvals)
    dur_avg, dur_rms = avg_and_rms(dur)
    
    #print('Npeaks:',Npeaks)
    #print('peaktimes:',peaktimes)
    #print('peakvals_avg:',peakvals_avg)
    #print('peakvals_rms:',peakvals_rms)
    #print('dur_avg:',dur_avg)
    #print('dur_rms:',dur_rms)
    return Npeaks, peaktimes, peakvals_avg,  peakvals_rms, dur_avg, dur_rms

if __name__ == '__main__':
    spikedurat = -40
    idur       = 1000 #100 # ms
    idelay     = 100
    iamp       = 0.4 # nA
    v_init     = -65 # mV
    Ra         = 100
    somasize   = 10 # 15 # 
    dendlen    = 1000
    denddiam   = 2
    nsegments  = 200 
    
    varymech = 'Na' # 'K' # 'leak'
    varyE_bool = True
    varyE = 50 #[30,40,50,60,70] #[30,40,70]# Every Cmf has 50 and 60. Need to run again for the other values
    varyg = 'None'
    
    varylist = [] # Should be redundant
    plotstring = '_vary'
    if varyE_bool==True:
        varylist = varyE
        plotstring = plotstring + 'E'
    else:
        varylist = varyg
        plotstring = plotstring + 'g'
      
    if varymech=='Na':
        folderstring = 'VaryNa/' 
        plotstring = plotstring + '_Na'
    elif varymech=='leak':
        folderstring = 'VaryLeak/'
        plotstring = plotstring + '_leak'
    elif varymech=='K':
        folderstring = 'VaryK/'
        plotstring = plotstring + '_K'

    changestring =''
    if varyE_bool==True:
        changestring =changestring+'_E'+str(varyE)+'_gdflt'
    else:
        changestring =changestring+'_Edf_g'+str(varyg)
    
    cm = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0]
    
    NCms = len(cm)
    
    Nspikes = numpy.zeros(NCms)
    avg_AP_ampl = numpy.zeros(NCms)
    rms_AP_ampl = numpy.zeros(NCms)
    avg_AP_halfwidth = numpy.zeros(NCms)
    rms_AP_halfwidth = numpy.zeros(NCms)
    
    # Set names
    outfolder = 'Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+folderstring
    currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    outfolder = outfolder+currentfolder
    outfilename_Nspikes = outfolder+'basHH_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmfs_Nspikes_vs_Cmsprx.txt'
    outfilename_APampl  = outfolder+'basHH_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmfs_Vmax_vs_Cmsprx.txt'
    outfilename_APdhw   = outfolder+'basHH_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmfs_sdurat%s_vs_Cmsprx.txt' % str(spikedurat)
    plotname_Nspikes    = outfolder+'basHH_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmfs_Nspikes_vs_Cmsprx.png'
    plotname_APampl     = outfolder+'basHH_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmfs_Vmax_vs_Cmsprx.png'
    plotname_APdhw      = outfolder+'basHH_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmfs_sdurat%s_vs_Cmsprx.png' % str(spikedurat)
    # make files
    outfile_Nspikes = open(outfilename_Nspikes,'w')
    outfile_APampl  = open(outfilename_APampl,'w')
    outfile_APdhw   = open(outfilename_APdhw,'w')
    for j in range(NCms):
        print('Step ', j+1, ' of', NCms)
        infolder = 'Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+ folderstring+currentfolder
        filename = infolder+'basHH_cmf'+str(cm[j])+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+changestring+'_V_sprxf.txt' 
        #try:
        #print('In try')
        Nspikes[j], peaktimes, avg_AP_ampl[j], rms_AP_ampl[j], avg_AP_halfwidth[j], rms_AP_halfwidth[j]= manual(filename,idelay,idur,spikedurat)
        #print('Nspikes[j]:',Nspikes[j])
        #if Nspikes[j]==1: # Poor statistics, omit this point
        #   continue
        #except:
        #    print('In except')
        #    continue
        outfile_Nspikes.write('%.5f %i\n' % (cm[j],Nspikes[j]))
        outfile_APampl.write('%.5f %.10f %.10f\n' % (cm[j],avg_AP_ampl[j],rms_AP_ampl[j]))
        outfile_APdhw.write('%.5f %.10f %.10f\n' % (cm[j],avg_AP_halfwidth[j],rms_AP_halfwidth[j]))
    outfile_Nspikes.close()
    outfile_APampl.close()
    outfile_APdhw.close()
    
    # Plot results
    plt.figure(figsize=(6,5))
    plt.plot(cm,Nspikes)
    plt.xlabel(r'$C_{m}$ [$\mu$ F/cm$^2$]')
    plt.ylabel(r'$N_{spikes}$')
    plt.title(r'Capacitance vs number of spikes')
    plt.tight_layout()
    plt.savefig(plotname_Nspikes)
    
    plt.figure(figsize=(6,5))
    plt.errorbar(cm,avg_AP_ampl, yerr=rms_AP_ampl, capsize=2)
    plt.xlabel(r'$C_{m}$ [$\mu$ F/cm$^2$]')
    plt.ylabel(r'Peak voltage [mV]')
    plt.title(r'Capacitance vs peak voltage')
    plt.tight_layout()
    plt.savefig(plotname_APampl)
    
    plt.figure(figsize=(6,5))
    plt.errorbar(cm,avg_AP_halfwidth, yerr=rms_AP_halfwidth, capsize=2)
    plt.xlabel(r'$C_{m} $[$\mu$ F/cm$^2$]')
    plt.ylabel(r'Spike duration at %s mV [ms] % str(spikedurat)')
    plt.title(r'Capacitance vs spike duration at %s mV' % str(spikedurat))
    plt.tight_layout()
    plt.savefig(plotname_APdhw)
    
    
    # Print results to terminal
    print('Nspikes:', Nspikes)
    print('AP amplitude, avg:', avg_AP_ampl)
    print('AP amplitude, rms:', rms_AP_ampl)
    print('AP duration at half width, avg:', avg_AP_halfwidth)
    print('AP duration at half width, rms:', rms_AP_halfwidth)