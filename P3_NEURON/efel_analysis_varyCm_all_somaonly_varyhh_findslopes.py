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
                                           ['peak_time','AP_amplitude', 'AP_duration_half_width', 'Spikecount', 'voltage_base'])
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
    # treat data and perform avg,rms where needed
    avg_AP_ampl, rms_AP_ampl = avg_and_rms(trace_results["AP_amplitude"])
    avg_AP_halfwidth, rms_AP_halfwidth = avg_and_rms(trace_results["AP_duration_half_width"])
    Nspikes = trace_results["Spikecount"]
    Nspikes = Nspikes[0]
    return Nspikes, avg_AP_ampl, rms_AP_ampl, avg_AP_halfwidth, rms_AP_halfwidth


if __name__ == '__main__':
    testmodel = 496497595 # 488462965 #
    idur      = 1000 # ms
    idelay    = 1
    iamp      = 0.02 # nA
    v_init    = -70 #-86.5 # mV
    Ra        = 100 # -150
    somasize  = 10
    
    # Default HH values:
    ena = 50
    ek = -77
    el_hh = -54.3
    gnabar_hh = 0.12
    gkbar_hh = 0.036
    gl_hh = 0.0003
    
    ### Change HH values here: ####
    #ena = 20
    #ek = -70
    #el_hh = -70
    #gnabar_hh = 0.14
    #gkbar_hh = 0.036*2
    #gl_hh = 0.000003
    
    iamps = [0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18]
    cms = [0.8,0.9,1.0,1.1,1.2,1.25,1.3,1.4,1.5]
        
    NCms = len(cms)
    Namps = len(iamps)
    slopes = numpy.zeros(Namps)
    
    outfolder   = 'Results/IStim/Soma%i/Vary_iamp/'%somasize
    outfilename = outfolder+'Cmslopes_vs_iamp.txt'
    plotname    = outfolder+'Cmslopes_vs_iamp.png'
    outfile     = open(outfilename,'w')
    for k in range(Namps):
        iamp = iamps[k]
        Nspikes = numpy.zeros(NCms)
        avg_AP_ampl = numpy.zeros(NCms)
        rms_AP_ampl = numpy.zeros(NCms)
        avg_AP_halfwidth = numpy.zeros(NCms)
        rms_AP_halfwidth = numpy.zeros(NCms)
        
        # Set names
        hhstring = '_ena'+str(ena)+'_ek'+str(ek)+'_el'+str(el_hh)+'_gnabar'+str(gnabar_hh)+'_gkbar'+str(gkbar_hh)+'_gl'+str(gl_hh)
        folder = 'Results/IStim/Soma%i/current_idur'%somasize+str(idur)+'_iamp'+str(iamp)+'/'
        # make files
        for j in range(NCms):
            print('Step ', j+1, ' of', NCms)
            cm = cms[j]
            filename = folder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+hhstring+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_V.txt' 
            Nspikes[j], avg_AP_ampl[j], rms_AP_ampl[j], avg_AP_halfwidth[j], rms_AP_halfwidth[j] = main(filename,idelay,idur)
    
        a = (Nspikes[2]-Nspikes[0])/(cms[2]-cms[0])
        print('Slope, iamp=',iamp, ': ', a) # Loop this?
        slopes[k] = a
        outfile.write('%.2f %.2e\n' % (iamp,a))
outfile.close()

plt.figure(figsize=(6,5))
plt.plot(iamps,slopes,'-o')
plt.xlabel(r'$I$ (nA)')
plt.ylabel(r'Slope (Hz cm$^2$/$\mu$F)')
plt.savefig(plotname)
plt.show()
