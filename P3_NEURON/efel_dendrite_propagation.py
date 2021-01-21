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
                                           ['time_to_first_spike'])#,'peak_time','AP_amplitude', 'AP_duration_half_width', 'Spikecount', 'voltage_base'])

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
    #avg_AP_ampl, rms_AP_ampl = avg_and_rms(trace_results["AP_amplitude"])
    #avg_AP_halfwidth, rms_AP_halfwidth = avg_and_rms(trace_results["AP_duration_half_width"])
    time_firstspike = trace_results["time_to_first_spike"]
    #Nspikes = trace_results["Spikecount"]
    #Nspikes = Nspikes[0]
    
    print('trace_results:',trace_results)
    print('traces_results:',traces_results)
    return time_firstspike #Nspikes, avg_AP_ampl, rms_AP_ampl, avg_AP_halfwidth, rms_AP_halfwidth


if __name__ == '__main__':
    testmodel = 496497595
    idur = 1 # ms
    iamp = 1.0 # nA
    idelay = 1
    v_init = -86.5
    
    dendnr = 6
    idxs = [0,6,25,39,58] # Not a great way to do this, but the moment I do this for more dendrites the read-in is harder to automate.
    Nidx = len(idxs)
    
    # Defaulting to original values:
    # DO NOT TOUCH THESE!
    # SET THEM BELOW INSTEAD!
    cm_soma = 1.14805
    cm_dend = 9.98231
    cm_axon = 3.00603
    
    # Changing values of membrane capacitance:
    #cm_soma = 0.01
    #cm_dend = 0.01
    #cm_axon = 10.0
    
    folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/dendritepropagation/'
    
    firstspiketimes = numpy.zeros(Nidx)
    for j in range(Nidx):
        filename = folder+'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) +'_vinit'+str(v_init)+ '_addedRa_V_dendsec%i_dend%i.txt' % (j,dendnr)
        #time_firstspike, Nspikes, avg_AP_ampl, rms_AP_ampl, avg_AP_halfwidth, rms_AP_halfwidth = main(filename,idelay,idur)
        firstspiketimes[j] = main(filename,idelay,idur)
    print('time_firstspike:', time_firstspike)