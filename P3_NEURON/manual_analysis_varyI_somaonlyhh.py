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
    peakmins  = []
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
        if voltage[i-1]>voltage[i] and voltage[i+1]>voltage[i] and voltage[i]<vthr:
            peakmins.append(voltage[i])
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
    for i in range(10,Ndur-1):
        dur.append(passtimes_down[i]-passtimes_up[i])
    for i in range(1,len(peaktimes)):
        isi.append(peaktimes[i]-peaktimes[i-1])
    isi      = isi[9:-1]
    peakvals = peakvals[9:-1]
    peakmins = peakmins[9:-1]
    
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
    peakmins_avg, peakmins_rms = avg_and_rms(peakmins)
    peakvals_avg, peakvals_rms = avg_and_rms(peakvals)
    dur_avg, dur_rms = avg_and_rms(dur)
    isi_avg, isi_rms = avg_and_rms(isi)
    
    #print('Npeaks:',Npeaks)
    #print('peaktimes:',peaktimes)
    #print('peakvals_avg:',peakvals_avg)
    #print('peakvals_rms:',peakvals_rms)
    #print('dur_avg:',dur_avg)
    #print('dur_rms:',dur_rms)
    return Npeaks, peaktimes, peakmins_avg, peakmins_rms, peakvals_avg,  peakvals_rms, dur_avg, dur_rms, isi_avg, isi_rms, isi

if __name__ == '__main__':
    spikedurat = -40
    idur       = 1000 #100 # ms
    idelay     = 10
    v_init     = -70 #-65 # mV
    Ra         = 100
    somasize   = 10 # 15 # 
    dendlen    = 1000
    denddiam   = 1
    nsegments  = 200 
    
    # Default HH values:
    ena = 50
    ek = -77
    el_hh = -54.3
    gnabar_hh = 0.12
    gkbar_hh = 0.036
    gl_hh = 0.0003
    
    # Changing stuff:
    ek = ek-20
    #el_hh = el_hh+10
    #gnabar_hh *= 1.0
    
    varymech = 'Na' # 'leak' # 'K' # 
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
    
    cm = 1.0
    iamps = [0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5]
    
    Namps = len(iamps)
    
    Nspikes = numpy.zeros(Namps)
    avg_ISI = numpy.zeros(Namps)
    rms_ISI = numpy.zeros(Namps)
    avg_AP_ampl = numpy.zeros(Namps)
    rms_AP_ampl = numpy.zeros(Namps)
    avg_AP_mins = numpy.zeros(Namps)
    rms_AP_mins = numpy.zeros(Namps)
    avg_AP_halfwidth = numpy.zeros(Namps)
    rms_AP_halfwidth = numpy.zeros(Namps)
    
    # Set names
    outfolder = 'Results/IStim/Soma%i/' % somasize
    hhstring = '_ena'+str(ena)+'_ek'+str(ek)+'_el'+str(el_hh)+'_gnabar'+str(gnabar_hh)+'_gkbar'+str(gkbar_hh)+'_gl'+str(gl_hh)
    outfilename_Nspikes = outfolder+'somaonlyHH_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+hhstring+'_Nspikes_vs_I.txt'
    outfilename_APampl  = outfolder+'somaonlyHH_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+hhstring+'_Vmax_vs_I.txt'
    outfilename_APmins  = outfolder+'somaonlyHH_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+hhstring+'_Vmin_vs_I.txt'
    outfilename_APdhw   = outfolder+'somaonlyHH_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+hhstring+'_sdurat%s_vs_I.txt' % str(spikedurat)
    outfilename_ISI     = outfolder+'somaonlyHH_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+hhstring+'_ISI_vs_I.txt'
    plotname_Nspikes    = outfolder+'somaonlyHH_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+hhstring+'_Nspikes_vs_Cmsprx.png'
    plotname_APampl     = outfolder+'somaonlyHH_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+hhstring+'_Vmax_vs_I.png'
    plotname_APmins     = outfolder+'somaonlyHH_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+hhstring+'_Vmin_vs_I.png'
    plotname_APdhw      = outfolder+'somaonlyHH_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+hhstring+'_sdurat%s_vs_I.png' % str(spikedurat)
    plotname_ISI        = outfolder+'somaonlyHH_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+hhstring+'_ISI_vs_I.png'
    # make files
    outfile_Nspikes = open(outfilename_Nspikes,'w')
    outfile_APampl  = open(outfilename_APampl,'w')
    outfile_APmins  = open(outfilename_APmins,'w')
    outfile_APdhw   = open(outfilename_APdhw,'w')
    outfile_ISI     = open(outfilename_ISI,'w')
    firedbefore = False
    for j in range(Namps):
        iamp = iamps[j]
        print('Step ', j+1, ' of', Namps, 'cm:', cm, 'iamp:', iamp)
        infolder = outfolder + '/current_idur%i_iamp'%idur+str(iamp)+'/'
        filename = infolder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+hhstring+'_Ra%i_vinit' %Ra+str(v_init)+'_V.txt' 
        #try:
        #print('In try')
        Nspikes[j], peaktimes, avg_AP_mins[j], rms_AP_mins[j], avg_AP_ampl[j], rms_AP_ampl[j], avg_AP_halfwidth[j], rms_AP_halfwidth[j], avg_ISI[j], rms_ISI[j], ISI = manual(filename,idelay,idur,spikedurat)
        if Nspikes[j]!=0 and firedbefore==False:
            firedbefore = True
        if Nspikes[j]!=0 or firedbefore==False:
            outfile_Nspikes.write('%.5f %i\n' % (iamp,Nspikes[j]))
        if Nspikes[j]!=0:
            outfile_APampl.write('%.5f %.10f %.10f\n' % (iamp,avg_AP_ampl[j],rms_AP_ampl[j]))
            outfile_APmins.write('%.5f %.10f %.10f\n' % (iamp,avg_AP_mins[j],rms_AP_mins[j]))
            outfile_APdhw.write('%.5f %.10f %.10f\n' % (iamp,avg_AP_halfwidth[j],rms_AP_halfwidth[j]))
            outfile_ISI.write('%.5f %.10f %.10f\n' % (iamp,avg_ISI[j],rms_ISI[j]))
            # Write all ISIs:
            outfilename_ISI_all = infolder+'basPV_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmf'+str(cm)+'_ISIall_vs_I.txt'
            outfile_ISI_all = open(outfilename_ISI_all,'w')
            for k in range(len(ISI)):
                outfile_ISI_all.write('%.10f ' % ISI[k])
            outfile_ISI_all.close()
    outfile_Nspikes.close()
    outfile_APampl.close()
    outfile_APmins.close()
    outfile_APdhw.close()
    outfile_ISI.close()
    
    # Plot results
    plt.figure(figsize=(6,5))
    plt.plot(iamps,Nspikes)
    plt.xlabel(r'$I$ [nA]')
    plt.ylabel(r'$N_{spikes}$')
    plt.title(r'Input current vs number of spikes, one comp. HH')
    plt.tight_layout()
    plt.savefig(plotname_Nspikes)
    
    plt.figure(figsize=(6,5))
    plt.errorbar(iamps,avg_AP_ampl, yerr=rms_AP_ampl, capsize=2)
    plt.xlabel(r'$I$ [nA]')
    plt.ylabel(r'Peak voltage [mV]')
    plt.title(r'Input current vs peak voltage, one comp. HH')
    plt.tight_layout()
    plt.savefig(plotname_APampl)
    
    plt.figure(figsize=(6,5))
    plt.errorbar(iamps,avg_AP_mins, yerr=rms_AP_mins, capsize=2)
    plt.xlabel(r'$I$ [nA]')
    plt.ylabel(r'Peak minima [mV]')
    plt.title(r'Input current vs min peak value at I=%.2f, one comp. HH' % iamp)
    plt.tight_layout()
    plt.savefig(plotname_APmins)
    
    plt.figure(figsize=(6,5))
    plt.errorbar(iamps,avg_AP_halfwidth, yerr=rms_AP_halfwidth, capsize=2)
    plt.xlabel(r'$I$ [nA]')
    plt.ylabel(r'Spike duration at %s mV [ms] % str(spikedurat)')
    plt.title(r'Input current vs spike duration at %s mV, one comp. HH' % str(spikedurat))
    plt.tight_layout()
    plt.savefig(plotname_APdhw)

    plt.figure(figsize=(6,5))
    plt.errorbar(iamps,avg_ISI, yerr=rms_ISI, capsize=2)
    plt.xlabel(r'$I$ [nA]')
    plt.ylabel(r'Interspike interval [ms]')
    plt.title(r'Input current vs interspike interval at I=%.2f, one comp. HH')
    plt.tight_layout()
    plt.savefig(plotname_ISI)
    
    
    # Print results to terminal
    print('Nspikes:', Nspikes)
    print('AP amplitude, avg:', avg_AP_ampl)
    print('AP amplitude, rms:', rms_AP_ampl)
    print('AP duration at half width, avg:', avg_AP_halfwidth)
    print('AP duration at half width, rms:', rms_AP_halfwidth)