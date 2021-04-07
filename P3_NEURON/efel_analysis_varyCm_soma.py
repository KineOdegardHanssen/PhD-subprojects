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
    smallCms = True
    testmodel = 478513407#488462965#480633479#478513437# #496497595
    idur   = 1000 # ms
    idelay = 100
    iamp   = 0.41 # 0.17# nA # Should be 0.17 for 478513407
    v_init = -86.5 # mV
    
    # Defaulting to original values:
    # DO NOT TOUCH THESE!
    # SET THEM BELOW INSTEAD!
    if testmodel==496497595:
        if smallCms==True:
            cm_soma = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        else:
            cm_soma = [0.01,0.1,0.5,1.14805,2.0,3.0]
        cm_dend = 9.98231
        cm_axon = 3.00603
        v_init = -86.5 # mV
    elif testmodel==488462965:
        if smallCms==True:
            cm_soma = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        else:
            cm_soma = [0.01,0.1,0.5,1.0,2.0,3.31732779736,4.0,5.0] # Have 10 too, but flatlines
        cm_dend = 3.31732779736
        cm_axon = 3.31732779736
        v_init = -86.5 # mV
    elif testmodel==478513407: # Should have I=0.17
        if smallCms==True:
            cm_soma = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        else:
            cm_soma = [0.01,0.1,0.5,1.0,2.0,3.0]
        cm_dend = 1.0
        cm_axon = 1.0
        v_init = -83.7
    elif testmodel==480633479:
        cm_soma = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.704866,0.8,0.9,1.0]
        cm_dend = 0.704866 # 0.704866118957
        cm_axon = 0.704866 # 0.704866118957
        v_init = -96.8
    elif testmodel==478513437:
        if smallCms==True:
            cm_soma = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        else:
            cm_soma = [0.01,0.1,0.5,1.0,2.0,2.34539964752,3.0,4.0]
        cm_dend = 2.34539964752
        cm_axon = 2.34539964752
        v_init = -86.8
    
    NCms = len(cm_soma)
    
    Nspikes = numpy.zeros(NCms)
    avg_AP_ampl = numpy.zeros(NCms)
    rms_AP_ampl = numpy.zeros(NCms)
    avg_AP_halfwidth = numpy.zeros(NCms)
    rms_AP_halfwidth = numpy.zeros(NCms)
    
    # Set names
    folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp'% (testmodel,idur)+str(iamp)+'/Varycm_soma/'
    if smallCms==True:
        outfilename_Nspikes = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_Nspikes_vs_Cmsoma_smallCm.txt'
        outfilename_APampl  = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APampl_vs_Cmsoma_smallCm.txt'
        outfilename_APdhw   = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APdurhalfwidth_vs_Cmsoma_smallCm.txt'
        plotname_Nspikes    = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_Nspikes_vs_Cmsoma_smallCm.png'
        plotname_APampl     = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APampl_vs_Cmsoma_smallCm.png'
        plotname_APdhw      = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APdurhalfwidth_vs_Cmsoma_smallCm.png'
    else:
        outfilename_Nspikes = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_Nspikes_vs_Cmsoma.txt'
        outfilename_APampl  = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APampl_vs_Cmsoma.txt'
        outfilename_APdhw   = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APdurhalfwidth_vs_Cmsoma.txt'
        plotname_Nspikes    = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_Nspikes_vs_Cmsoma.png'
        plotname_APampl     = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APampl_vs_Cmsoma.png'
        plotname_APdhw      = folder+'cellmodel%i_current_idur%i_iamp'% (testmodel,idur)+str(iamp) +'_APdurhalfwidth_vs_Cmsoma.png'
    # make files
    outfile_Nspikes = open(outfilename_Nspikes,'w')
    outfile_APampl  = open(outfilename_APampl,'w')
    outfile_APdhw   = open(outfilename_APdhw,'w')
    for j in range(NCms):
        print('Step ', j+1, ' of', NCms)
        filename = folder +'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma[j]) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) + '_vinit'+str(v_init)+'_addedRa.txt'
        Nspikes[j], avg_AP_ampl[j], rms_AP_ampl[j], avg_AP_halfwidth[j], rms_AP_halfwidth[j] = main(filename,idelay,idur)
        outfile_Nspikes.write('%.5f %i\n' % (cm_soma[j],Nspikes[j]))
        outfile_APampl.write('%.5f %.10f %.10f\n' % (cm_soma[j],avg_AP_ampl[j],rms_AP_ampl[j]))
        outfile_APdhw.write('%.5f %.10f %.10f\n' % (cm_soma[j],avg_AP_halfwidth[j],rms_AP_halfwidth[j]))
    outfile_Nspikes.close()
    outfile_APampl.close()
    outfile_APdhw.close()
    # Plot results
    plt.figure(figsize=(6,5))
    plt.plot(cm_soma,Nspikes)
    plt.xlabel(r'$C_{m}$ of soma [$\mu$ F/cm$^2$]')
    plt.ylabel(r'$N_{spikes}$')
    plt.title(r'Soma capacitance vs number of spikes')
    plt.tight_layout()
    plt.savefig(plotname_Nspikes)
    
    plt.figure(figsize=(6,5))
    plt.errorbar(cm_soma,avg_AP_ampl, yerr=rms_AP_ampl, capsize=2)
    plt.xlabel(r'$C_{m}$ of soma [$\mu$ F/cm$^2$]')
    plt.ylabel(r'Spike amplitude [mV]')
    plt.title(r'Soma capacitance vs AP amplitude')
    plt.tight_layout()
    plt.savefig(plotname_APampl)
    
    plt.figure(figsize=(6,5))
    plt.errorbar(cm_soma,avg_AP_halfwidth, yerr=rms_AP_halfwidth, capsize=2)
    plt.xlabel(r'$C_{m}$ of soma [$\mu$ F/cm$^2$]')
    plt.ylabel(r'AP duration at half width [ms]')
    plt.title(r'Soma capacitance vs AP duration at half with')
    plt.tight_layout()
    plt.savefig(plotname_APdhw)
    
    
    
    
    # Print results to terminal
    print('Nspikes:', Nspikes)
    print('AP amplitude, avg:', avg_AP_ampl)
    print('AP amplitude, rms:', rms_AP_ampl)
    print('AP duration at half width, avg:', avg_AP_halfwidth)
    print('AP duration at half width, rms:', rms_AP_halfwidth)