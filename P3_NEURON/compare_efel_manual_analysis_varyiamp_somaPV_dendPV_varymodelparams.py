"""Modified from 'Basic example 1 for eFEL'."""

import efel
import numpy
import matplotlib.pyplot as plt

def avg_and_rms(x):
    try:
        N = len(x)
    except:
        avgx = 0
        rmsx = 0 
        return avgx,rmsx
    avgx = numpy.mean(x)
    rmsx = 0
    for i in range(N):
        rmsx += (avgx-x[i])**2
    rmsx = numpy.sqrt(rmsx/(N-1))
    return avgx,rmsx

def manual(filename,idelay,idur):
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
    vthr  = -20#vmax-0.15*deltav # If there is a peak above this value, we count it
    vprev = vthr-40 # A peak never kicks in at initiation, even if I change vthr
    durthr = -20 # Height at which we measure the duration. Need to be calibrated by peaks?
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
    
    #print(':',Npeaks)
    #print(':',peaktimes)
    #print(':',peakvals_avg)
    #print(':',peakvals_rms)
    #print(':',dur_avg)
    #print(':',dur_rms)
    return Npeaks, peaktimes, peakvals_avg,  peakvals_rms, dur_avg, dur_rms
    
    
    
    
def main(filename,idelay,idur):
    """Main"""

    # Use numpy to read the trace data from the txt file
    data = numpy.loadtxt(filename)

    # Time is the first column
    time = data[:, 0]
    # Voltage is the second column
    voltage = data[:, 1]
    
    # Now we will construct the datastructure that will be passed to eFEL

    # A 'trace' is a dictionary
    trace1 = {}

    # Set the 'T' (=time) key of the trace
    trace1['T'] = time

    # Set the 'V' (=voltage) key of the trace
    trace1['V'] = voltage

    # Set the 'stim_start' (time at which a stimulus starts, in ms)
    # key of the trace
    # Warning: this need to be a list (with one element)
    trace1['stim_start'] = [idelay]

    # Set the 'stim_end' (time at which a stimulus end) key of the trace
    # Warning: this need to be a list (with one element)
    trace1['stim_end'] = [idelay+idur]

    # Multiple traces can be passed to the eFEL at the same time, so the
    # argument should be a list
    traces = [trace1]

    # Now we pass 'traces' to the efel and ask it to calculate the feature
    # values
    traces_results = efel.getFeatureValues(traces,
                                           ['peak_time','peak_voltage','AP_amplitude', 'AP_duration_half_width', 'Spikecount', 'voltage_base'])
    ###### This is only printing. I do not really need it when I'm looping ##############
    '''
    # The return value is a list of trace_results, every trace_results
    # corresponds to one trace in the 'traces' list above (in same order)
    for trace_results in traces_results:
        # trace_result is a dictionary, with as keys the requested features
        for feature_name, feature_values in trace_results.items():
            #print('feature_values:',feature_values)
            if len(feature_values)!=0: # I changed this from if feature_values!=None:
                print("Feature %s has the following values: %s" % \
                (feature_name, ', '.join([str(x) for x in feature_values])))
    '''
    #####################################################################################
    trace_results = traces_results[0] # Because I am only looping over one cell, I guess
    APampl = trace_results["AP_amplitude"]
    APdurhw = trace_results["AP_duration_half_width"]
    Nspikes = trace_results["Spikecount"]
    Vpeak   = trace_results["peak_voltage"]
    print('APampl:',APampl)
    print('APdurhw:',APdurhw)
    print('Nspikes:',Nspikes)

    #------------------------ Basic data analysis --------------------------------------------------#
    # treat data and perform avg,rms where needed
    avg_AP_ampl, rms_AP_ampl = avg_and_rms(trace_results["AP_amplitude"])
    avg_AP_halfwidth, rms_AP_halfwidth = avg_and_rms(trace_results["AP_duration_half_width"])
    avg_Vpeak, rms_Vpeak = avg_and_rms(trace_results["peak_voltage"])
    Nspikes = trace_results["Spikecount"]
    Nspikes = Nspikes[0]
    Nampl   = Nspikes
    #--------------------------- Safeguards --------------------------------------------------------#
    ## Duration of APs: Two conditions need to be met: AP dur under 3ms and no more than 3stdv>avg
    AP_dur_ok = []
    basic_AP_dur = trace_results["AP_duration_half_width"]
    rms_threshold = avg_AP_halfwidth+3.*rms_AP_halfwidth
    try:
        for i in range(len(basic_AP_dur)):
            tempdur = basic_AP_dur[i]
            if tempdur<3 and tempdur<rms_threshold: # Not too large deviation from the rest
                AP_dur_ok.append(tempdur)
    except:
        AP_dur_ok = []
    AP_dur_ok = numpy.array(AP_dur_ok)
    avg_AP_halfwidth, rms_AP_halfwidth = avg_and_rms(AP_dur_ok)
    Ndur = len(AP_dur_ok)
    ## First spike gets overestimated amplitude (happens to every sim., but best not to use bad data)
    ampl_data = trace_results["AP_amplitude"]
    try:
        if len(ampl_data)>1: # Don't want to throw away everything
            avg_AP_ampl, rms_AP_ampl = avg_and_rms(ampl_data[1:]) # Skipping the first peak
            Nampl -=1    
    except:
        Nampl=Nampl
    return Nspikes, avg_AP_ampl, rms_AP_ampl, avg_AP_halfwidth, rms_AP_halfwidth, Nampl, Ndur, avg_Vpeak, rms_Vpeak


if __name__ == '__main__':
    # varymech: Difficult to make this make sense for other ions than Na (and leak)...
    # 478513407: gbar_Kv2like=0.0854285
    varymech = 'NaV' #'Im_v2'#'Kv2like'#'Kd'#'pas' #
    varyE = -60 # Default...407Kd-107 Default478513407pas: -83.6528;  Default478513407Na: 53
    varyg = 0.9 # Default...407Kd2.94396e-010 # Default478513407pas: 0.000362109; Default478513407Na: 0.0409177
    if varymech=='NaV':
        folderstring = 'VaryNa/' 
        textsnippet = '_varyNa'
    elif varymech=='pas':
        folderstring = 'VaryPas/'
        textsnippet = '_varypas'
    elif varymech=='Kd':
        folderstring = 'VaryKd/'
        textsnippet = '_varyKd'
    elif varymech=='Kv2like':
        folderstring = 'VaryKv2like/'
        textsnippet = '_varyKv2like'
    elif varymech=='Kv3_1':
        folderstring = 'VaryKv3_1/'
        textsnippet = '_varyKv3_1'
    elif varymech=='SK':
        folderstring = 'VarySK/'
        textsnippet = '_varySK'
    elif varymech=='K_T':
        folderstring = 'VaryK_T/'
        textsnippet = '_varyK_T'
    elif varymech=='Im_v2':
        folderstring = 'VaryIm_v2/'
        textsnippet = '_varyIm_v2'
    #
    testmodel = 478513407#488462965#478513437#
    idur   = 1000 #100 # ms
    idelay = 100
    iamps  = [0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18] # nA
    v_init = -86.5 # mV
    Ra      = 150

    cm = 1.0
    
    N = len(iamps)
    
    Nspikes = numpy.zeros(N)
    avg_AP_ampl = numpy.zeros(N)
    rms_AP_ampl = numpy.zeros(N)
    avg_AP_halfwidth = numpy.zeros(N)
    rms_AP_halfwidth = numpy.zeros(N)
    avg_Vpeak = numpy.zeros(N)
    rms_Vpeak = numpy.zeros(N)
    # Manually 
    Npeaks_man    = numpy.zeros(N)
    #peaktimes_man = numpy.zeros(N) 
    peakvals_man_avg  = numpy.zeros(N)
    dur_man_avg       = numpy.zeros(N)
    peakvals_man_rms  = numpy.zeros(N)
    dur_man_rms       = numpy.zeros(N)
    
    # Set names
    folder_base = 'Results/%i/IStim/' % testmodel + folderstring
    outfilename_Nspikes = folder_base+'basps_cellmodel%i_current_idur%i'% (testmodel,idur)+textsnippet+'_E'+str(varyE)+'_g'+str(varyg)+'_Nspikes_vs_iamp_compare.txt'
    plotname_Nspikes    = folder_base+'basps_cellmodel%i_current_idur%i'% (testmodel,idur)+textsnippet+'_E'+str(varyE)+'_g'+str(varyg)+'_Nspikes_vs_iamp_compare.png'
    plotname_APampl     = folder_base+'basps_cellmodel%i_current_idur%i'% (testmodel,idur)+textsnippet+'_E'+str(varyE)+'_g'+str(varyg)+'_Vpeak_vs_iamp_compare.png'
    plotname_APdhw      = folder_base+'basps_cellmodel%i_current_idur%i'% (testmodel,idur)+textsnippet+'_E'+str(varyE)+'_g'+str(varyg)+'_APdurhalfwidth_vs_iamp_compare.png'
    outfile = open(outfilename_Nspikes,'w')
    for j in range(N):
        iamp = iamps[j]
        print('Step ', j+1, ' of', N)
        currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
        folder = folder_base+currentfolder
        filename = folder+'basps_cm'+str(cm)+'_E'+str(varyE)+'_g'+str(varyg)+'_Ra%i_vinit' %Ra+str(v_init)+'_V.txt' 
        Nspikes[j], avg_AP_ampl[j], rms_AP_ampl[j], avg_AP_halfwidth[j], rms_AP_halfwidth[j], Nampl, Ndur, avg_Vpeak[j], rms_Vpeak[j] = main(filename,idelay,idur)
        #Npeaks_man[j], peaktimes_man, peakvals_man_avg[j], peakvals_man_rms[j], dur_man_avg[j], dur_man_rms[j] = manual(filename,idelay,idur)
        outfile.write('%.2f %i\n' % (iamp,Nspikes[j]))
    outfile.close()
    
    # Plot results
    plt.figure(figsize=(6,5))
    plt.plot(iamps,Nspikes, '-o', label='eFEL')
    #plt.plot(iamps,Npeaks_man, '-o', label='Manually')
    plt.xlabel(r'$I$ [nA]')
    plt.ylabel(r'$N_{spikes}$')
    plt.title(r'Capacitance vs number of spikes')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(plotname_Nspikes)
    
    plt.figure(figsize=(6,5))
    plt.errorbar(iamps,avg_Vpeak, yerr=rms_Vpeak, capsize=2, label='eFEL')
    #plt.errorbar(iamps,peakvals_man_avg, yerr=peakvals_man_rms, capsize=2, label='Manually')
    plt.xlabel(r'$I$ [nA]')
    plt.ylabel(r'Peak voltage [mV]')
    plt.title(r'Capacitance vs peak voltage')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(plotname_APampl)
    
    plt.figure(figsize=(6,5))
    plt.errorbar(iamps,avg_AP_halfwidth, yerr=rms_AP_halfwidth, capsize=2, label='eFEL')
    #plt.errorbar(iamps,dur_man_avg, yerr=dur_man_rms, capsize=2, label='Manually')
    plt.xlabel(r'$I$ [nA]')
    plt.ylabel(r'AP width at half amplitude [ms]')
    plt.title(r'Capacitance vs AP width')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(plotname_APdhw)
    
    plt.show()
    
    # Print results to terminal
    print('Nspikes:', Nspikes)
    print('AP amplitude, avg:', avg_AP_ampl)
    print('AP amplitude, rms:', rms_AP_ampl)
    print('AP duration at half width, avg:', avg_AP_halfwidth)
    print('AP duration at half width, rms:', rms_AP_halfwidth)