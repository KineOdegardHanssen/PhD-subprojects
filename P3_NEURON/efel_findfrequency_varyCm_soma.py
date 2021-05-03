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
                                           ['peak_time', 'AP_amplitude', 'AP_duration_half_width',  'Spikecount', 'voltage_base'])
    # Keep some for testing
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
    trace_results = traces_results[0] # Because I am only looping over one cell
    #------------------------ Basic data analysis --------------------------------------------------#
    # treat data and perform avg,rms where needed
    Nspikes = trace_results["Spikecount"]
    Nspikes = Nspikes[0]
    
    # For double checking: 
    AP_dur_data  = trace_results["AP_duration_half_width"]
    ampl_data    = trace_results["AP_amplitude"]
    peak_time    = trace_results["peak_time"]
    voltage_base = trace_results["voltage_base"]
    
    '''
    print('Nspikes:', Nspikes, '; len(peak_time):',len(peak_time), '; len(AP_dur_data):' , len(AP_dur_data), '; ampl_data:', len(ampl_data))
    print('voltage_base:',voltage_base)
    print('peak_time:', peak_time)
    
    print('AP_dur:', AP_dur_data)
    print('ampl:', ampl_data)
    
    plt.figure(figsize=(6,5))
    plt.plot(time,voltage)
    plt.xlabel('Time [ms]')
    plt.ylabel('Voltage [mV]')
    plt.legend(loc='upper right')
    plt.show()
    '''
          
    return Nspikes


if __name__ == '__main__':
    smallCms = True#False#
    testmodel = 478513407#478513437#488462965# # Bad:#480633479# #496497595
    idur   = 1000 # ms
    idelay = 100
    iamps  = [0,0.02,0.04,0.06,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,
0.38,0.4,0.42,0.44,0.46,0.48,0.5]
    v_init = -86.5 # mV
    idur_SI = idur/1000. # Time is measured in ms
    Ni = len(iamps)
    
    # Defaulting to original values:
    # DO NOT TOUCH THESE!
    # SET THEM BELOW INSTEAD!
    if testmodel==496497595:
        cm_soma = 1.14805
        cm_dend = 9.98231
        cm_axon = 3.00603
        v_init = -86.5 # mV
    elif testmodel==488462965:
        cm_soma = 3.31732779736
        cm_dend = 3.31732779736
        cm_axon = 3.31732779736
        v_init = -86.5 # mV
    elif testmodel==478513407: # Should have I=0.17
        cm_soma = 1.0
        cm_dend = 1.0
        cm_axon = 1.0
        v_init = -83.7
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
    cm_soma = 0.5  
    
    fs = numpy.zeros(Ni)
    
    # Set names
    folder = 'Allen_test_changecapacitance/figures/%i/' % testmodel
    if smallCms==True:
        outfilename = folder+'f_I_curves/f_I_curves_cmsoma'+str(cm_soma)+'.txt'
        plotname    = folder+'f_I_curves/f_I_curves_cmsoma'+str(cm_soma)+'.png'
    # make files
    outfile = open(outfilename,'w')
    for j in range(Ni):
        iamp = iamps[j]
        print('Step ', j+1, ' of', Ni)
        infolder = folder+'current_idur%i_iamp'% idur + str(iamp) + '/Varycm_soma/'
        filename = infolder +'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) + '_vinit'+str(v_init)+'_addedRa.txt'
        try:
            Nspikes = main(filename,idelay,idur)
        except:
            infolder = folder+'current_idur%i_iamp'% idur + str(iamp) + '/'
            filename = infolder +'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) + '_vinit'+str(v_init)+'_addedRa.txt'
            Nspikes = main(filename,idelay,idur)
        f = Nspikes/float(idur_SI)
        fs[j] = f
        outfile.write('%.5f %.5f\n' % (iamp,f))
        #print('f:',f)
    outfile.close()
    # Plot results
    plt.figure(figsize=(6,5))
    plt.plot(iamps,fs,'-o')
    plt.xlabel(r'Input current [nA]')
    plt.ylabel(r'Frequence [Hz]')
    plt.title(r'Soma capacitance vs number of spikes, Cm soma %.2f $\mu$F/cm$^2$' % cm_soma)
    plt.tight_layout()
    plt.savefig(plotname)
    