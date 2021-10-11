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
    avg_AP_halfwidth, rms_AP_halfwidth = avg_and_rms(AP_dur_ok)
    Ndur = len(AP_dur_ok)
    ## First spike gets overestimated amplitude (happens to every sim., but best not to use bad data)
    ampl_data = trace_results["AP_amplitude"]
    if len(ampl_data)>1: # Don't want to throw away everything
        avg_AP_ampl, rms_AP_ampl = avg_and_rms(ampl_data[1:]) # Skipping the first peak
        Nampl -=1    
    return Nspikes, avg_AP_ampl, rms_AP_ampl, avg_AP_halfwidth, rms_AP_halfwidth, Nampl, Ndur


if __name__ == '__main__':
    idur      = 1000 #100 # ms
    idelay    = 100
    iamp      = 0.4 # nA
    v_init    = -65 # mV
    Ra        = 100
    somasize  = 10 # 15 # 
    dendlen   = 1000
    denddiam  = 2
    nsegments = 200 
    
    varymech = 'Na' # 'K' # 'leak'
    varyE_bool = True
    varyE = 50 #[30,40,50,60,70] #[30,40,70]# Every Cmf has 50 and 60. Need to run again for the other values
    varyg = 'None'
    
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
    elif varymech=='leak':
        folderstring = 'VaryLeak/'
        plotstring = plotstring + '_leak'
    elif varymech=='K':
        folderstring = 'VaryK/'
        plotstring = plotstring + '_K'

    changestring =''
    if varyE_bool==True:
        changestring =changestring+'_E'+str(varyE)+'_gdflt'
    else:
        changestring =changestring+'_Edf_g'+str(varyg)
    
    cm = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0]
    
    NCms = len(cm)
    
    Nspikes = numpy.zeros(NCms)
    avg_AP_ampl = numpy.zeros(NCms)
    rms_AP_ampl = numpy.zeros(NCms)
    avg_AP_halfwidth = numpy.zeros(NCms)
    rms_AP_halfwidth = numpy.zeros(NCms)
    
    # Set names
    outfolder = 'Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+folderstring
    currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    outfolder = outfolder+currentfolder
    outfilename_Nspikes = outfolder+'basHH_current_idur%i_iamp'% (idur)+str(iamp) +'_cmfs_Nspikes_vs_Cmsprx.txt'
    outfilename_APampl  = outfolder+'basHH_current_idur%i_iamp'% (idur)+str(iamp) +'_cmfs_APampl_vs_Cmsprx.txt'
    outfilename_APdhw   = outfolder+'basHH_current_idur%i_iamp'% (idur)+str(iamp) +'_cmfs_APdurhalfwidth_vs_Cmsprx.txt'
    plotname_Nspikes    = outfolder+'basHH_current_idur%i_iamp'% (idur)+str(iamp) +'_cmfs_Nspikes_vs_Cmsprx.png'
    plotname_APampl     = outfolder+'basHH_current_idur%i_iamp'% (idur)+str(iamp) +'_cmfs_APampl_vs_Cmsprx.png'
    plotname_APdhw      = outfolder+'basHH_current_idur%i_iamp'% (idur)+str(iamp) +'_cmfs_APdurhalfwidth_vs_Cmsprx.png'
    # make files
    outfile_Nspikes = open(outfilename_Nspikes,'w')
    outfile_APampl  = open(outfilename_APampl,'w')
    outfile_APdhw   = open(outfilename_APdhw,'w')
    for j in range(NCms):
        print('Step ', j+1, ' of', NCms)
        infolder = 'Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+ folderstring+currentfolder
        filename = infolder+'basHH_cmf'+str(cm[j])+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+changestring+'_V_sprxf.txt' 
        try:
            print('In try')
            Nspikes[j], avg_AP_ampl[j], rms_AP_ampl[j], avg_AP_halfwidth[j], rms_AP_halfwidth[j], Nampl, Ndur = main(filename,idelay,idur)
            #print('Nspikes[j]:',Nspikes[j])
            #if Nspikes[j]==1: # Poor statistics, omit this point
            #   continue
        except:
            print('In except')
            continue
        outfile_Nspikes.write('%.5f %i\n' % (cm[j],Nspikes[j]))
        outfile_APampl.write('%.5f %.10f %.10f\n' % (cm[j],avg_AP_ampl[j],rms_AP_ampl[j]))
        outfile_APdhw.write('%.5f %.10f %.10f\n' % (cm[j],avg_AP_halfwidth[j],rms_AP_halfwidth[j]))
    outfile_Nspikes.close()
    outfile_APampl.close()
    outfile_APdhw.close()
    
    # Plot results
    plt.figure(figsize=(6,5))
    plt.plot(cm,Nspikes)
    plt.xlabel(r'$C_{m}$ [$\mu$ F/cm$^2$]')
    plt.ylabel(r'$N_{spikes}$')
    plt.title(r'Capacitance vs number of spikes')
    plt.tight_layout()
    plt.savefig(plotname_Nspikes)
    
    plt.figure(figsize=(6,5))
    plt.errorbar(cm,avg_AP_ampl, yerr=rms_AP_ampl, capsize=2)
    plt.xlabel(r'$C_{m}$ [$\mu$ F/cm$^2$]')
    plt.ylabel(r'Spike amplitude [mV]')
    plt.title(r'Capacitance vs AP amplitude')
    plt.tight_layout()
    plt.savefig(plotname_APampl)
    
    plt.figure(figsize=(6,5))
    plt.errorbar(cm,avg_AP_halfwidth, yerr=rms_AP_halfwidth, capsize=2)
    plt.xlabel(r'$C_{m} $[$\mu$ F/cm$^2$]')
    plt.ylabel(r'AP width at half amplitude [ms]')
    plt.title(r'Capacitance vs AP width at half amplitude')
    plt.tight_layout()
    plt.savefig(plotname_APdhw)
    
    
    # Print results to terminal
    print('Nspikes:', Nspikes)
    print('AP amplitude, avg:', avg_AP_ampl)
    print('AP amplitude, rms:', rms_AP_ampl)
    print('AP duration at half width, avg:', avg_AP_halfwidth)
    print('AP duration at half width, rms:', rms_AP_halfwidth)