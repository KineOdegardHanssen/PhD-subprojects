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
        if time[i]>idur+idelay:
            break
    # Find resting potential:
    k = 0
    
    print('Nspikes before testing:',Npeaks)
    # Checking if we've got consistent firing:
    if Npeaks!=0:
        if peaktimes[-1]<=(idur/2+idelay): # Checking if there's no firing in the last half of the stim. interval
            Npeaks=0                       # Is that a proper limit? Last third instead?
    if Npeaks==0:
        passtimes_up   = []
        passtimes_down = []
    
    # I should probably plot stuff, to make sure... If-test?
    dur = []
    isi = []
    Ndur = min([len(passtimes_up),len(passtimes_down)])
    for i in range(10,Ndur-1):
        dur.append(passtimes_down[i]-passtimes_up[i])
    for i in range(1,len(peaktimes)):
        thisisi = peaktimes[i]-peaktimes[i-1]
        isi.append(thisisi)
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
    
    #print('dur:',dur)
    
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
    cm         = 1.0
    spikedurat = -40
    idur       = 1000 #100 # ms
    idelay     = 100
    v_init     = -65 # mV
    Ra         = 100
    somasize   = 10 # 15 # 
    dendlen    = 1000
    denddiam   = 1
    nsegments  = 200 
    gna        = 0.12
    gfactor    = 2.5 # [1.5,2.0,3.0,5.0]#[0.8,1.0,1.2]#[0.1,0.5]#
    
    varymech = 'Na' # 'K' # 'leak'
    varyE_bool = False # True # 
    varyg_bool = True # False # 
    varyE = 50 #[30,40,50,60,70] #[30,40,70]# Every Cmf has 50 and 60. Need to run again for the other values
    varyg = gfactor*gna
    
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
        namestring = 'na'
    elif varymech=='leak':
        folderstring = 'VaryLeak/'
        plotstring = plotstring + '_leak'
        namestring = 'l'
    elif varymech=='K':
        folderstring = 'VaryK/'
        plotstring = plotstring + '_K'
        namestring = 'k'

    changestring =''
    if varyE_bool==True:
        changestring =changestring+'_E'+str(varyE)+'_gdf'
    else:
        changestring =changestring+'_g'+str(varyg)
    
    iamps     = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5]
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
    outfolder = 'Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+folderstring
    outfilename_Nspikes = outfolder+'basHHdpas_manual_cm'+str(cm)+'_g'+namestring+str(varyg)+'_Nspikes_vs_Iamp.txt'
    outfilename_APampl  = outfolder+'basHHdpas_manual_cm'+str(cm)+'_g'+namestring+str(varyg)+'_Vmax_vs_Iamp.txt'
    outfilename_APmins  = outfolder+'basHHdpas_manual_cm'+str(cm)+'_g'+namestring+str(varyg)+'_Vmin_vs_Iamp.txt'
    outfilename_APdhw   = outfolder+'basHHdpas_manual_cm'+str(cm)+'_g'+namestring+'s_sdurat%s_vs_Iamp.txt' % str(spikedurat)
    outfilename_ISI     = outfolder+'basHHdpas_manual_cm'+str(cm)+'_g'+namestring+str(varyg)+'_ISI_vs_Iamp.txt'
    plotname_Nspikes    = outfolder+'basHHdpas_manual_cm'+str(cm)+'_g'+namestring+str(varyg)+'_Nspikes_vs_Iamp.png'
    plotname_APampl     = outfolder+'basHHdpas_manual_cm'+str(cm)+'_gs_Vmax_vs_Iamp.png'
    plotname_APmins     = outfolder+'basHHdpas_manual_cm'+str(cm)+'_g'+namestring+str(varyg)+'_Vmin_vs_Iamp.png'
    plotname_APdhw      = outfolder+'basHHdpa_manual_cm'+str(cm)+'_g'+namestring+'s_sdurat%s_vs_Iamp.png' % str(spikedurat)
    plotname_ISI        = outfolder+'basHHdpas_manual_cm'+str(cm)+'_g'+namestring+str(varyg)+'_ISI_vs_Iamp.png'
    # make files
    outfile_Nspikes = open(outfilename_Nspikes,'w')
    outfile_APampl  = open(outfilename_APampl,'w')
    outfile_APmins  = open(outfilename_APmins,'w')
    outfile_APdhw   = open(outfilename_APdhw,'w')
    outfile_ISI     = open(outfilename_ISI,'w')
    firedbefore = False
    for j in range(Namps):
        iamp = iamps[j]
        print('Step ', j+1, ' of', Namps, ', iamp=',iamp)
        currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
        infolder = 'Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+ folderstring+currentfolder
        filename = infolder+'basHHdpas_cmf'+str(cm)+'_Ra%i_vinit' %Ra+str(v_init)+changestring+'_V_all.txt' 
        #try:
        #print('In try')
        Nspikes[j], peaktimes, avg_AP_mins[j], rms_AP_mins[j], avg_AP_ampl[j], rms_AP_ampl[j], avg_AP_halfwidth[j], rms_AP_halfwidth[j], avg_ISI[j], rms_ISI[j], ISI = manual(filename,idelay,idur,spikedurat)
        #print('Nspikes[j]:',Nspikes[j])
        #if Nspikes[j]==1: # Poor statistics, omit this point
        #   continue
        #except:
        #    print('In except')
        #    continue
        if Nspikes[j]!=0 and firedbefore==False:
            firedbefore = True
        if Nspikes[j]!=0 or firedbefore==False:
            outfile_Nspikes.write('%.2f %i\n' % (iamp,Nspikes[j]))
        if Nspikes[j]!=0:
            outfile_APampl.write('%.2f %.10f %.10f\n' % (iamp,avg_AP_ampl[j],rms_AP_ampl[j]))
            outfile_APmins.write('%.2f %.10f %.10f\n' % (iamp,avg_AP_mins[j],rms_AP_mins[j]))
            outfile_APdhw.write('%.2f %.10f %.10f\n' % (iamp,avg_AP_halfwidth[j],rms_AP_halfwidth[j]))
            outfile_ISI.write('%.5f %.10f %.10f\n' % (iamp,avg_ISI[j],rms_ISI[j]))
            # Write all ISIs:
            outfilename_ISI_all = outfolder+'basHHdpas_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmf'+str(cm)+'_g'+namestring+str(varyg)+'_ISIall_vs_Iamp.txt'
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
    plt.xlabel(r'$C_{m}$ [$\mu$ F/cm$^2$]')
    plt.ylabel(r'$N_{spikes}$')
    plt.title(r'Capacitance vs number of spikes, BASHH')
    plt.tight_layout()
    plt.savefig(plotname_Nspikes)
    
    plt.figure(figsize=(6,5))
    plt.errorbar(iamps,avg_AP_ampl, yerr=rms_AP_ampl, capsize=2)
    plt.xlabel(r'$C_{m}$ [$\mu$ F/cm$^2$]')
    plt.ylabel(r'Peak voltage [mV]')
    plt.title(r'Capacitance vs peak voltage, BASHH')
    plt.tight_layout()
    plt.savefig(plotname_APampl)
    
    plt.figure(figsize=(6,5))
    plt.errorbar(iamps,avg_AP_mins, yerr=rms_AP_mins, capsize=2)
    plt.xlabel(r'$C_{m}$ [$\mu$ F/cm$^2$]')
    plt.ylabel(r'Peak voltage [mV]')
    plt.title(r'Capacitance vs min peak value at I=%.2f, BAS HH' % iamp)
    plt.tight_layout()
    plt.savefig(plotname_APmins)
    
    plt.figure(figsize=(6,5))
    plt.errorbar(iamps,avg_AP_halfwidth, yerr=rms_AP_halfwidth, capsize=2)
    plt.xlabel(r'$C_{m} $[$\mu$ F/cm$^2$]')
    plt.ylabel(r'Spike duration at %s mV [ms] % str(spikedurat)')
    plt.title(r'Capacitance vs spike duration at %s mV, BASHH' % str(spikedurat))
    plt.tight_layout()
    plt.savefig(plotname_APdhw)

    plt.figure(figsize=(6,5))
    plt.errorbar(iamps,avg_ISI, yerr=rms_ISI, capsize=2)
    plt.xlabel(r'$C_{m} $[$\mu$ F/cm$^2$]')
    plt.ylabel(r'Interspike interval [ms]')
    plt.title(r'Capacitance vs interspike interval at I=%.2f, BASHH' % iamp)
    plt.tight_layout()
    plt.savefig(plotname_ISI)
    
    
    # Print results to terminal
    #print('Nspikes:', Nspikes)
    #print('AP amplitude, avg:', avg_AP_ampl)
    #print('AP amplitude, rms:', rms_AP_ampl)
    #print('AP duration at half width, avg:', avg_AP_halfwidth)
    #print('AP duration at half width, rms:', rms_AP_halfwidth)