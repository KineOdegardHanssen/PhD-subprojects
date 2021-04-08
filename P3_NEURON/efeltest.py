"""Modified from 'Basic example 1 for eFEL'."""

import efel
import numpy

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
    
    print('len(time):', len(time))
    print('len(voltage):', len(voltage))
    print('time:', time)
    print('voltage:', voltage)
    
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
                                           ['peak_time','AP_amplitude', 'AP_duration_half_width', 'Spikecount', 'voltage_base'])

    # The return value is a list of trace_results, every trace_results
    # corresponds to one trace in the 'traces' list above (in same order)
    for trace_results in traces_results:
        # trace_result is a dictionary, with as keys the requested features
        for feature_name, feature_values in trace_results.items():
            print('feature_values:',feature_values)
            if len(feature_values)!=0: # I changed this from if feature_values!=None:
                print("Feature %s has the following values: %s" % \
                (feature_name, ', '.join([str(x) for x in feature_values])))
    # treat data and perform avg,rms where needed
    # I think I will have separate script for plotting V of different Cm's together.
    print('trace_results:',trace_results)
    #print('trace_results["AP_amplitude"]):',trace_results["AP_amplitude"])
    avg_AP_ampl, rms_AP_ampl = avg_and_rms(trace_results["AP_amplitude"])
    avg_AP_halfwidth, rms_AP_halfwidth = avg_and_rms(trace_results["AP_duration_half_width"])
    Nspikes = trace_results["Spikecount"]
    Nspikes = Nspikes[0]
    #--------------------------- Safeguards --------------------------------------------------------#
    # Duration of APs: Two conditions need to be met: AP dur under 3ms and no more than 3stdv>avg
    AP_dur_ok = []
    basic_AP_dur = trace_results["AP_duration_half_width"]
    rms_threshold = avg_AP_halfwidth+3.*rms_AP_halfwidth
    for i in range(len(basic_AP_dur)):
        tempdur = basic_AP_dur[i]
        if tempdur<3. and tempdur<rms_threshold: # Not too large deviation from the rest
            print('tempdur:',tempdur, '; rms threshold:',rms_threshold)
            AP_dur_ok.append(tempdur)
    AP_dur_ok = numpy.array(AP_dur_ok)
    avg_AP_ampl, rms_AP_ampl = avg_and_rms(AP_dur_ok)
    
    print('trace_results:',trace_results)
    #print('traces_results:',traces_results)
    print('AP_dur_ok:',AP_dur_ok)
    print('avg_AP_halfwidth:',avg_AP_halfwidth, '; rms_AP_halfwidth:',rms_AP_halfwidth)
    return Nspikes, avg_AP_ampl, rms_AP_ampl, avg_AP_halfwidth, rms_AP_halfwidth


if __name__ == '__main__':
    testmodel = 480633479#478513407#478513437#488462965# #496497595
    idur = 1000  #2  # ms
    iamp = 0.41 #1.0 #0.5 # nA
    idelay = 100 #1  #
    v_init = -86.5
    
    # Defaulting to original values:
    # DO NOT TOUCH THESE!
    # SET THEM BELOW INSTEAD!
    if testmodel==496497595:
        cm_soma = 1.14805
        cm_dend = 9.98231
        cm_axon = 3.00603
    elif testmodel==488462965:
        cm_soma = 3.31732779736
        cm_dend = 3.31732779736
        cm_axon = 3.31732779736
    elif testmodel==480633479:
        cm_soma = 0.704866
        cm_dend = 0.704866 # 0.704866118957
        cm_axon = 0.704866 # 0.704866118957
        v_init = -96.8
    elif testmodel==478513437:
        cm_soma = 2.34539964752
        cm_dend = 2.34539964752
        cm_axon = 2.34539964752
        v_init = -86.8
    elif testmodel==478513407:
        cm_soma = 1.0
        cm_dend = 1.0
        cm_axon = 1.0
        v_init = -83.7
    
    
    # Changing values of membrane capacitance:
    #cm_soma = 0.9#0.4
    #cm_dend = cm_soma# 5.0
    #cm_axon = cm_soma# 5.0
    
    
    folder = 'Allen_test_changecapacitance/'
    filename = folder+'figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) + '_vinit'+str(v_init)+'_addedRa.txt'
    Nspikes, avg_AP_ampl, rms_AP_ampl, avg_AP_halfwidth, rms_AP_halfwidth = main(filename,idelay,idur)
    print('Nspikes:', Nspikes)
    print('AP amplitude:', avg_AP_ampl, ' +/-', rms_AP_ampl)
    print('AP duration at half width:', avg_AP_halfwidth, ' +/-', rms_AP_halfwidth)