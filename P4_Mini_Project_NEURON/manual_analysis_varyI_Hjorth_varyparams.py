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
    passtimes_up = []
    passvals_up  = []
    passtimes_down = []
    passvals_down  = []
    for i in range (1,len(voltage)-1):  
        if time[i]<idur+idelay:
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
        else:
            break
    
    # Checking if we've got consistent firing:
    if Npeaks!=0:
        if peaktimes[-1]<=(idur/2.+idelay): #Checking if there's no firing in the last half of the stim. interval
            Npeaks=0                        # Is that a proper limit? Last third instead?
    
    # I should probably plot stuff, to make sure... If-test?
    dur = []
    isi = []
    Ndur = min([len(passtimes_up),len(passtimes_down)])
    for i in range(2,Ndur-1):
        dur.append(passtimes_down[i]-passtimes_up[i])
    for i in range(1,len(peaktimes)):
        isi.append(peaktimes[i]-peaktimes[i-1])
    isi      = isi[1:-1]
    peakvals = peakvals[1:-1]
    peakmins = peakmins[1:-1]
    time_peakvals = peaktimes[1:-1]
    
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
    idelay     = 100
    v_init     = -86 # mV
    Ra         = 100
    somasize   = 10 # 15 # 
    dtexp      = -7
    t_before_rec = -600.
    
    
    model_folder = '' # We are in the model folder
    
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
    
    vary_naf     = True # False
    vary_kaf     = True # False
    vary_SK      = False # 
    vary_Ca_HVA  = False
    vary_gpas    = False # 
    
    gnaf   = 0.5
    gkaf   = 0.5
    gsk    = 1.0
    gcahva = 1.0
    gpas   = 1.0
    
    if vary_naf==True:
        namestring = namestring + '_gnaf'+str(gnav)+'p'
    if vary_kaf==True:
        namestring = namestring + '_gkaf'+str(gkaf)+'p'
    if vary_SK==True:
        namestring = namestring + '_gSK'+str(gsk)+'p'
    if vary_Ca_HVA==True:
        namestring = namestring + '_gCaHVA'+str(gcahva)+'p'
    if vary_gpas==True: 
        namestring = namestring + '_gpas'+str(gpas)+'p'
    namestring = namestring +'_'
    
    
    cm = 1.0
    iamps = [0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03,0.031,0.032,0.033,0.034,0.035,0.036,0.037,0.038,0.039,0.04,0.041,0.042,0.043,0.044,0.045,0.046,0.047,0.048,0.049,0.05,0.051,0.052,0.053,0.054,0.055,0.056,0.057,0.058,0.059,0.06]    
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
    outfolder = 'Results/Soma%i/' % somasize
    outfilename_Nspikes = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Nspikes_vs_I.txt'
    outfilename_APampl  = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Vmax_vs_I.txt'
    outfilename_APmins  = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Vmin_vs_I.txt'
    outfilename_APdhw   = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_sdurat%s_vs_I.txt' % str(spikedurat)
    outfilename_ISI     = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_ISI_vs_I.txt'
    plotname_Nspikes    = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Nspikes_vs_I.png'
    plotname_APampl     = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Vmax_vs_I.png'
    plotname_APmins     = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Vmin_vs_I.png'
    plotname_APdhw      = outfolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_sdurat%s_vs_I.png' % str(spikedurat)
    # make files
    outfile_Nspikes = open(outfilename_Nspikes,'w')
    outfile_APampl  = open(outfilename_APampl,'w')
    outfile_APmins  = open(outfilename_APmins,'w')
    outfile_APdhw   = open(outfilename_APdhw,'w')
    firedbefore = False
    for j in range(Namps):
        iamp = iamps[j]
        print('Step ', j+1, ' of', Namps, 'cm:', cm, 'iamp:', iamp)
        infolder = outfolder + 'current_idur%i_iamp'%idur+str(iamp)+'/'
        filename = infolder+namestring+'somaonly_cm'+str(cm_factor)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_V.txt'
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
            # Write all durs:
            outfilename_durs_all = infolder+'somaHjorth_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmf'+str(cm)+namestring+str(v_init)+'_trec'+str(t_before_rec)+'_alldurs.txt'
            figname_durs_all     = infolder+'somaHjorth_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmf'+str(cm)+namestring+str(v_init)+'_trec'+str(t_before_rec)+'_alldurs.png'
            outfile_durs_all = open(outfilename_durs_all,'w')
            for k in range(len(dur)):
                outfile_durs_all.write('%.10f ' % dur[k])
            outfile_durs_all.close()
            ptimes = peaktimes[1:-1]
            if len(ptimes)==len(dur):
                plt.figure(figsize=(6,5))
                plt.plot(ptimes,dur)
                plt.xlabel('Time [ms]')
                plt.ylabel('Voltage [mV]')
                plt.title('I=%s nA' % str(iamp))
                plt.tight_layout()
                plt.savefig(figname_durs_all)
            # Write all peak times (matching durs):
            outfilename_pt_all = infolder+'somaHjorth_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmf'+str(cm)+namestring+str(v_init)+'_trec'+str(t_before_rec)+'_peaktimes.txt'
            outfile_pt_all = open(outfilename_pt_all,'w')
            for k in range(len(peaktimes)-1):
                outfile_pt_all.write('%.10f ' % peaktimes[k])
            outfile_pt_all.close()
            #outfile_ISI.write('%.5f %.10f %.10f\n' % (iamp,avg_ISI[j],rms_ISI[j]))
            # Write all ISIs:
            #outfilename_ISI_all = infolder+'basPV_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmf'+str(cm)+'_ISIall_vs_I.txt'
            #outfile_ISI_all = open(outfilename_ISI_all,'w')
            #for k in range(len(ISI)):
            #    outfile_ISI_all.write('%.10f ' % ISI[k])
            #outfile_ISI_all.close()
    outfile_Nspikes.close()
    outfile_APampl.close()
    outfile_APmins.close()
    outfile_APdhw.close()
    #outfile_ISI.close()
    
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
    
    '''
    plt.figure(figsize=(6,5))
    plt.errorbar(iamps,avg_ISI, yerr=rms_ISI, capsize=2)
    plt.xlabel(r'$I$ [nA]')
    plt.ylabel(r'Interspike interval [ms]')
    plt.title(r'Input current vs interspike interval at I=%.2f, one comp. HH')
    plt.tight_layout()
    plt.savefig(plotname_ISI)
    '''
    
    
    # Print results to terminal
    print('Nspikes:', Nspikes)
    print('AP amplitude, avg:', avg_AP_ampl)
    print('AP amplitude, rms:', rms_AP_ampl)
    print('AP duration at half width, avg:', avg_AP_halfwidth)
    print('AP duration at half width, rms:', rms_AP_halfwidth)
