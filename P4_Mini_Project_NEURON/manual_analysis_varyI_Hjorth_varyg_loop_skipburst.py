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

def manual(filename,idelay,idur,spikedurat,skiptime):
    """Manual approach"""

    # Use numpy to read the trace data from the txt file
    data = numpy.loadtxt(filename)

    # Time is the first column
    time = data[:, 0]
    # Voltage is the second column
    voltage = data[:, 1]
    
    anfraction = (idur-skiptime)/float(idur)
    tstartan   = idelay+skiptime
    print('anfraction:',anfraction)
    print('tstartan:',tstartan)
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
    passtimes_up = []
    passvals_up  = []
    passtimes_down = []
    passvals_down  = []
    minalready = False
    for i in range (1,len(voltage)-1):  
        if time[i]<idelay+idur:
            if voltage[i-1]<voltage[i] and voltage[i+1]<voltage[i] and voltage[i]>vthr and time[i]>tstartan:
                peaktimes.append(time[i])
                peakvals.append(voltage[i])
                Npeaks+=1
                minalready = False
            if voltage[i-1]>voltage[i] and voltage[i+1]>voltage[i] and voltage[i]<vthr and time[i]>tstartan and minalready==False and Npeaks>0:
                peakmins.append(voltage[i])
                minalready = True
            if voltage[i]>=durthr and voltage[i-1]<durthr and time[i]>tstartan: # Passing upwards
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
            elif voltage[i]>=durthr and voltage[i+1]<durthr and len(passtimes_up)>0: # Passing downwards
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
        else:
            break
    
    Npeaks /= anfraction # Want frequency
    # Checking if we've got consistent firing:
    if Npeaks!=0:
        if peaktimes[-1]<=(idur/2.+idelay): #Checking if there's no firing in the last half of the stim. interval
            Npeaks=0                        # Is that a proper limit? Last third instead?
    
    # I should probably plot stuff, to make sure... If-test?
    dur = []
    isi = []
    amps = []
    Namps = min([len(peakmins),len(peakvals)])
    Ndur = min([len(passtimes_up),len(passtimes_down)]) # Should be the same
    for i in range(Ndur-1):
        dur.append(passtimes_down[i]-passtimes_up[i])
    for i in range(Namps):
        amps.append(peakvals[i]-peakmins[i])
    for i in range(1,len(peaktimes)):
        isi.append(peaktimes[i]-peaktimes[i-1])
    #isi      = isi[:-1]
    #peakvals = peakvals[:-1]
    #peakmins = peakmins[:-1]
    time_peakvals = peaktimes#[:-1]
    
    '''
    plt.figure(figsize=(6,5))
    plt.plot(time,voltage,',')
    plt.plot(time_peakvals,peakvals,'o',label='peaks')
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
    amps_avg, amps_rms = avg_and_rms(amps)
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
    return Npeaks, peaktimes, peakmins_avg, peakmins_rms, peakvals_avg,  peakvals_rms, dur_avg, dur_rms, isi_avg, isi_rms, isi, dur, amps_avg, amps_rms

if __name__ == '__main__':
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
    
    gnaf    = 1.0
    gkaf    = 1.0
    gpas    = 1.0
    
    cm = 1.0
    iamps =  [0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03,0.031,0.032,0.033,0.034,0.035,0.036,0.037,0.038,0.039,0.04,0.041,0.042,0.043,0.044,0.045,0.046,0.047,0.048,0.049,0.05,0.051,0.052,0.053,0.054,0.055,0.056,0.057,0.058,0.059,0.06]    
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
    avg_AP_amplitudes = numpy.zeros(Namps)
    rms_AP_amplitudes = numpy.zeros(Namps)
    
    # Set names
    outfolder = 'Results/Soma%i/' % somasize
    outfilename_Nspikes = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'
    outfilename_amps    = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Ampl_vs_I_s'+str(skiptime)+'.txt'
    outfilename_APampl  = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Vmax_vs_I_s'+str(skiptime)+'.txt'
    outfilename_APmins  = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Vmin_vs_I_s'+str(skiptime)+'.txt'
    outfilename_APdhw   = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_sdurat%s' % str(spikedurat)+'_vs_I_s'+str(skiptime)+'.txt'
    outfilename_ISI     = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_ISI_vs_I_s'+str(skiptime)+'.txt'
    # make files
    outfile_Nspikes = open(outfilename_Nspikes,'w')
    outfile_amps    = open(outfilename_amps,'w')
    outfile_APampl  = open(outfilename_APampl,'w')
    outfile_APmins  = open(outfilename_APmins,'w')
    outfile_APdhw   = open(outfilename_APdhw,'w')
    outfile_ISI     = open(outfilename_ISI,'w')
    firedbefore = False
    for j in range(Namps):
        iamp = iamps[j]
        print('Step ', j+1, ' of', Namps, 'cm:', cm, 'iamp:', iamp)
        infolder = outfolder + 'current_idur%i_iamp'%idur+str(iamp)+'/'
        filename = infolder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_V.txt'
        Nspikes[j], peaktimes, avg_AP_mins[j], rms_AP_mins[j], avg_AP_ampl[j], rms_AP_ampl[j], avg_AP_halfwidth[j], rms_AP_halfwidth[j], avg_ISI[j], rms_ISI[j], ISI, durs, avg_AP_amplitudes[j], rms_AP_amplitudes[j] = manual(filename,idelay,idur,spikedurat,skiptime)
        if Nspikes[j]!=0 and firedbefore==False:
            firedbefore = True    
        if Nspikes[j]!=0 or firedbefore==False:
            outfile_Nspikes.write('%.5f %i\n' % (iamp,Nspikes[j]))
        if Nspikes[j]!=0:
            outfile_amps.write('%.5f %.10f %.10f\n' % (iamp,avg_AP_amplitudes[j],rms_AP_amplitudes[j]))
            outfile_APampl.write('%.5f %.10f %.10f\n' % (iamp,avg_AP_ampl[j],rms_AP_ampl[j]))
            outfile_APmins.write('%.5f %.10f %.10f\n' % (iamp,avg_AP_mins[j],rms_AP_mins[j]))
            outfile_APdhw.write('%.5f %.10f %.10f\n' % (iamp,avg_AP_halfwidth[j],rms_AP_halfwidth[j]))
            # Write all durs:
            outfilename_durs_all = infolder+'somaHjorth_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmf'+str(cm)+str(v_init)+'_trec'+str(t_before_rec)+'_alldurs_s'+str(skiptime)+'.txt'
            figname_durs_all     = infolder+'somaHjorth_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmf'+str(cm)+str(v_init)+'_trec'+str(t_before_rec)+'_alldurs'+str(skiptime)+'.png'
            outfile_durs_all = open(outfilename_durs_all,'w')
            for k in range(len(durs)):    
                outfile_durs_all.write('%.10f ' % durs[k])
            outfile_durs_all.close()
            ptimes = peaktimes[1:-1]
            if len(ptimes)==len(durs):
                plt.figure(figsize=(6,5))
                plt.plot(ptimes,durs)
                plt.xlabel('Time [ms]')
                plt.ylabel('Voltage [mV]')    
                plt.title('I=%s nA' % str(iamp))
                plt.tight_layout()
                plt.savefig(figname_durs_all)
            # Write all peak times (matching durs):    
            outfilename_pt_all = infolder+'somaHjorth_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmf'+str(cm)+str(v_init)+'_trec'+str(t_before_rec)+'_peaktimes'+str(skiptime)+'.txt'
            outfile_pt_all = open(outfilename_pt_all,'w')
            for k in range(len(peaktimes)-1):    
                outfile_pt_all.write('%.10f ' % peaktimes[k])
            outfile_pt_all.close()
            outfile_ISI.write('%.5f %.10f %.10f\n' % (iamp,avg_ISI[j],rms_ISI[j]))
            # Write all ISIs:
            outfilename_ISI_all = infolder+'somaHjorth_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmf'+str(cm)+str(v_init)+'_trec'+str(t_before_rec)+'_ISIall'+str(skiptime)+'.txt'
            outfile_ISI_all = open(outfilename_ISI_all,'w')
            for k in range(len(ISI)):
                outfile_ISI_all.write('%.10f ' % ISI[k])
            outfile_ISI_all.close()
    outfile_Nspikes.close()
    outfile_amps.close()
    outfile_APampl.close()
    outfile_APmins.close()
    outfile_APdhw.close()
    outfile_ISI.close()

