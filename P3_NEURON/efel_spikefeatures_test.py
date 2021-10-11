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

def main(filename,idelay,idur):
    """Main"""

    # Use numpy to read the trace data from the txt file
    data = numpy.loadtxt(filename)

    # Time is the first column
    time = data[:, 0]
    # Voltage is the second column
    voltage = data[:, 1]
    
    '''
    plt.figure(figsize=(6,5))
    plt.plot(time,voltage)
    plt.xlabel('t [ms]')
    plt.ylabel('V [mV]')
    plt.title('V vs t')
    plt.tight_layout()
    plt.show()
    '''
    
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
                                           ['peak_time','AP_amplitude', 'AP_duration_half_width', 'Spikecount', 'voltage_base','AP_begin_indices', 'peak_voltage'])
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
    #------------------------ Basic data analysis --------------------------------------------------#
    # treat data and perform avg,rms where needed
    avg_AP_ampl, rms_AP_ampl = avg_and_rms(trace_results["AP_amplitude"])
    avg_AP_halfwidth, rms_AP_halfwidth = avg_and_rms(trace_results["AP_duration_half_width"])
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
    print('basic_AP_dur:',basic_AP_dur)
    print('AP_dur_ok:',AP_dur_ok)
    avg_AP_halfwidth, rms_AP_halfwidth = avg_and_rms(AP_dur_ok)
    Ndur = len(AP_dur_ok)
    ## First spike gets overestimated amplitude (happens to every sim., but best not to use bad data)
    ampl_data = trace_results["AP_amplitude"]
    if len(ampl_data)>1: # Don't want to throw away everything
        avg_AP_ampl, rms_AP_ampl = avg_and_rms(ampl_data[1:]) # Skipping the first peak
        Nampl -=1    
    print('ampl_data[1:]:',ampl_data[1:])
    print('ampl_data:',ampl_data)
    #
    
    beg_in = trace_results["AP_begin_indices"]
    print('Voltage at start:',voltage[beg_in])
    print('Average voltage at start:',numpy.mean(voltage[beg_in]))
    print('Max voltage at start:',max(voltage[beg_in]))
    print('Min voltage at start:',min(voltage[beg_in]))
    pv = trace_results["peak_voltage"]
    print('Peaks:',pv)
    print('avg(Peaks):',numpy.mean(pv))
    vgbi = voltage[beg_in]
    #vgbi = vgbi[1:]
    print('avg(Peaks-Voltage at start):',numpy.mean(pv-vgbi))
    print('avg(Peaks)-avg(Voltage at start)):',numpy.mean(pv)-numpy.mean(voltage[beg_in]))
    print('Peaks-Voltage at start-ampl_data:',pv-vgbi-ampl_data)
    print('numpy.mean(Peaks-Voltage at start-ampl_data):',numpy.mean(pv-vgbi-ampl_data))
    return Nspikes, avg_AP_ampl, rms_AP_ampl, avg_AP_halfwidth, rms_AP_halfwidth, Nampl, Ndur


if __name__ == '__main__':
    testmodel = 478513437#488462965 #496497595 # 
    idur   = 1000 #100 # ms
    idelay = 100
    iamp   = 0.5#0.41 # nA
    v_init = -86.5 # mV
    Ra      = 150
    
    # Defaulting to original values:
    # DO NOT TOUCH THESE!
    # SET THEM BELOW INSTEAD!
    cm = []
    if testmodel==496497595:
        cm = [0.5,1,2,3,4,4.5,5,5.5,6,7]#,8,9,10]
    elif testmodel==488462965:
        cm = [0.8]#[0.5,1,2,3,4,4.5,5,5.5,6,7]#,8,9,10] # Not made
    elif testmodel==478513437:
        cm = [0.7]
    cm = [0.8]
    NCms = len(cm)
    
    Nspikes = numpy.zeros(NCms)
    avg_AP_ampl = numpy.zeros(NCms)
    rms_AP_ampl = numpy.zeros(NCms)
    avg_AP_halfwidth = numpy.zeros(NCms)
    rms_AP_halfwidth = numpy.zeros(NCms)
    
    # Set names
    outfolder = 'Results/%i/IStim/' % testmodel
    currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    outfolder = outfolder+currentfolder
    for j in range(NCms):
        print('Step ', j+1, ' of', NCms)
        filename = outfolder+'basps_cm'+str(cm[j])+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_V.txt' 
        Nspikes[j], avg_AP_ampl[j], rms_AP_ampl[j], avg_AP_halfwidth[j], rms_AP_halfwidth[j], Nampl, Ndur = main(filename,idelay,idur)
    
    '''
    # Plot results
    plt.figure(figsize=(6,5))
    plt.plot(cm,Nspikes)
    plt.xlabel(r'$C_{m}$ [$\mu$ F/cm$^2$]')
    plt.ylabel(r'$N_{spikes}$')
    plt.title(r'Capacitance vs number of spikes')
    plt.tight_layout()
    plt.show()
    
    plt.figure(figsize=(6,5))
    plt.errorbar(cm,avg_AP_ampl, yerr=rms_AP_ampl, capsize=2)
    plt.xlabel(r'$C_{m}$ [$\mu$ F/cm$^2$]')
    plt.ylabel(r'Spike amplitude [mV]')
    plt.title(r'Capacitance vs AP amplitude')
    plt.tight_layout()
    plt.show()
    
    plt.figure(figsize=(6,5))
    plt.errorbar(cm,avg_AP_halfwidth, yerr=rms_AP_halfwidth, capsize=2)
    plt.xlabel(r'$C_{m} $[$\mu$ F/cm$^2$]')
    plt.ylabel(r'AP width at half amplitude [ms]')
    plt.title(r'Capacitance vs AP width at half amplitude')
    plt.tight_layout()
    plt.show()
    '''
    
    
    # Print results to terminal
    print('Nspikes:', Nspikes)
    print('AP amplitude, avg:', avg_AP_ampl)
    print('AP amplitude, rms:', rms_AP_ampl)
    print('AP duration at half width, avg:', avg_AP_halfwidth)
    print('AP duration at half width, rms:', rms_AP_halfwidth)
    print('cm:',cm)