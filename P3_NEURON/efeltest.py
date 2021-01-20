"""Modified from 'Basic example 1 for eFEL'."""

import efel
import numpy

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
                                           ['AP_begin_indices','AP_amplitude', 'AP_duration_half_width', 'voltage_base'])

    # The return value is a list of trace_results, every trace_results
    # corresponds to one trace in the 'traces' list above (in same order)
    for trace_results in traces_results:
        # trace_result is a dictionary, with as keys the requested features
        for feature_name, feature_values in trace_results.items():
            print('feature_values:',feature_values)
            if len(feature_values)!=0: # I changed this from if feature_values!=None:
                print("Feature %s has the following values: %s" % \
                (feature_name, ', '.join([str(x) for x in feature_values])))


if __name__ == '__main__':
    idelay = 1
    idur   = 100
    
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
    
    testmodel = 496497595
    idur = 100 # ms
    iamp = 0.5 # nA
    folder = 'Allen_test_changecapacitance/'
    filename = folder+'figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) + '.txt'
    main(filename,idelay,idur)