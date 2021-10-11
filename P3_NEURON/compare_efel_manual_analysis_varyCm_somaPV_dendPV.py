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
    for i in range(len(basic_AP_dur)):
        tempdur = basic_AP_dur[i]
        if tempdur<3 and tempdur<rms_threshold: # Not too large deviation from the rest
            AP_dur_ok.append(tempdur)
    AP_dur_ok = numpy.array(AP_dur_ok)
    avg_AP_halfwidth, rms_AP_halfwidth = avg_and_rms(AP_dur_ok)
    Ndur = len(AP_dur_ok)
    ## First spike gets overestimated amplitude (happens to every sim., but best not to use bad data)
    ampl_data = trace_results["AP_amplitude"]
    if len(ampl_data)>1: # Don't want to throw away everything
        avg_AP_ampl, rms_AP_ampl = avg_and_rms(ampl_data[1:]) # Skipping the first peak
        Nampl -=1    
    return Nspikes, avg_AP_ampl, rms_AP_ampl, avg_AP_halfwidth, rms_AP_halfwidth, Nampl, Ndur, avg_Vpeak, rms_Vpeak


if __name__ == '__main__':
    testmodel = 478513407#478513437#488462965#
    idur   = 1000 #100 # ms
    idelay = 100
    iamp   = 0.3 # nA
    v_init = -86.5 # mV
    Ra      = 150
    
    # Defaulting to original values:
    # DO NOT TOUCH THESE!
    # SET THEM BELOW INSTEAD!
    cm = [0.5, 1, 2, 10]
    # Have smallCm instead?:
    cm = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0] #[0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0]#
    
    NCms = len(cm)
    
    Nspikes = numpy.zeros(NCms)
    avg_AP_ampl = numpy.zeros(NCms)
    rms_AP_ampl = numpy.zeros(NCms)
    avg_AP_halfwidth = numpy.zeros(NCms)
    rms_AP_halfwidth = numpy.zeros(NCms)
    avg_Vpeak = numpy.zeros(NCms)
    rms_Vpeak = numpy.zeros(NCms)
    # Manually 
    Npeaks_man    = numpy.zeros(NCms)
    #peaktimes_man = numpy.zeros(NCms) 
    peakvals_man_avg  = numpy.zeros(NCms)
    dur_man_avg       = numpy.zeros(NCms)
    peakvals_man_rms  = numpy.zeros(NCms)
    dur_man_rms       = numpy.zeros(NCms)
    
    # Set names
    outfolder = 'Results/%i/IStim/' % testmodel
    currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    outfolder = outfolder+currentfolder
    plotname_Nspikes    = outfolder+'basps_cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_Nspikes_vs_Cmall_compare.png'
    plotname_APampl     = outfolder+'basps_cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_Vpeak_vs_Cmall_compare.png'
    plotname_APdhw      = outfolder+'basps_cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APdurhalfwidth_vs_Cmall_compare.png'
    for j in range(NCms):
        print('Step ', j+1, ' of', NCms)
        filename = outfolder+'basps_cm'+str(cm[j])+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_V.txt' 
        Nspikes[j], avg_AP_ampl[j], rms_AP_ampl[j], avg_AP_halfwidth[j], rms_AP_halfwidth[j], Nampl, Ndur, avg_Vpeak[j], rms_Vpeak[j] = main(filename,idelay,idur)
        Npeaks_man[j], peaktimes_man, peakvals_man_avg[j], peakvals_man_rms[j], dur_man_avg[j], dur_man_rms[j] = manual(filename,idelay,idur)
    
    # Plot results
    plt.figure(figsize=(6,5))
    plt.plot(cm,Nspikes, '-o', label='eFEL')
    plt.plot(cm,Npeaks_man, '-o', label='Manually')
    plt.xlabel(r'$C_{m}$ [$\mu$ F/cm$^2$]')
    plt.ylabel(r'$N_{spikes}$')
    plt.title(r'Capacitance vs number of spikes')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(plotname_Nspikes)
    
    plt.figure(figsize=(6,5))
    plt.errorbar(cm,avg_Vpeak, yerr=rms_Vpeak, capsize=2, label='eFEL')
    plt.errorbar(cm,peakvals_man_avg, yerr=peakvals_man_rms, capsize=2, label='Manually')
    plt.xlabel(r'$C_{m}$ [$\mu$ F/cm$^2$]')
    plt.ylabel(r'Peak voltage [mV]')
    plt.title(r'Capacitance vs peak voltage')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(plotname_APampl)
    
    plt.figure(figsize=(6,5))
    plt.errorbar(cm,avg_AP_halfwidth, yerr=rms_AP_halfwidth, capsize=2, label='eFEL')
    plt.errorbar(cm,dur_man_avg, yerr=dur_man_rms, capsize=2, label='Manually')
    plt.xlabel(r'$C_{m} $[$\mu$ F/cm$^2$]')
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