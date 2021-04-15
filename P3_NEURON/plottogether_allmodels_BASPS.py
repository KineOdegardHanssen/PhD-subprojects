"""Modified from 'Basic example 1 for eFEL'."""

import efel
import numpy
import numpy as np
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
    trace_results = traces_results[0] # Because I am only looping over one cell
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
    idur   = 1000 # ms
    idelay = 100
    iamp   = 0.6 # nA
    Ra     = 150
    v_init = -86.5 #-65 # mV
    
    smallCms = True
    if smallCms==True:
        cms = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0] #[0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.6,0.8,0.9,1.0]
    else:
        cms = [0.5,1,2,3,4,4.5,5,5.5,6,7,8,9,10]
    testmodels = [488462965,478513407,478513437] #[488462965,478513407,480633479,478513437]
    NCms = len(cms)
    Ncms = NCms  # Typo a lot of places.
    Nmodels = len(testmodels)
    
    outfolder      = 'Results/Comparemodels/'
    avgandrmsname  = outfolder+'compiled_avgandrms_BAPS.png'
    everymodelname = outfolder+'compiled_alltogether_BAPS.png'
    
    # Model-specific arrays:
    Nspikes_bymodel      = []
    avg_AP_ampl_bymodel      = []
    rms_AP_ampl_bymodel      = []
    avg_AP_halfwidth_bymodel = []
    rms_AP_halfwidth_bymodel = []
    # I will do averaging in the same go
    Nspikes_avg          = np.zeros(Ncms)
    Nampl_avg            = np.zeros(Ncms)
    Ndur_avg             = np.zeros(Ncms)
    AP_ampl_avg          = np.zeros(Ncms)
    AP_halfwidth_avg     = np.zeros(Ncms)
    Nspikes_rms          = np.zeros(Ncms)
    AP_ampl_rms          = np.zeros(Ncms)
    AP_halfwidth_rms     = np.zeros(Ncms)
    Nspikes_all          = np.zeros((Ncms,Nmodels)) # j,i
    Nampl_all            = np.zeros((Ncms,Nmodels)) # j,i
    Ndur_all             = np.zeros((Ncms,Nmodels)) # j,i
    AP_ampl_all          = np.zeros((Ncms,Nmodels))
    AP_halfwidth_all     = np.zeros((Ncms,Nmodels))
    AP_ampl_rms_all      = np.zeros((Ncms,Nmodels))
    AP_halfwidth_rms_all = np.zeros((Ncms,Nmodels))
    
    folder = 'Results/Comparemodels/'
    plotname_Nspikes    = folder+'basps_idur%i_iamp'% (idur)+str(iamp) +'_Nspikes_vs_Cm'
    plotname_APampl     = folder+'basps_idur%i_iamp'% (idur)+str(iamp) +'_APampl_vs_Cm'
    plotname_APdhw      = folder+'basps_idur%i_iamp'% (idur)+str(iamp) +'_APdurhalfwidth_vs_Cm'
    plotname_Nspikes_avgrms = folder+'basps_idur%i_iamp'% (idur)+str(iamp) +'_Nspikes_vs_Cm'
    plotname_APampl_avgrms  = folder+'basps_idur%i_iamp'% (idur)+str(iamp) +'_APampl_vs_Cm'
    plotname_APdhw_avgrms   = folder+'basps_idur%i_iamp'% (idur)+str(iamp) +'_APdurhalfwidth_vs_Cm'
    if smallCms==True:
        plotname_Nspikes    = plotname_Nspikes +'_smallCm.png'
        plotname_APampl     = plotname_APampl +'_smallCm.png'
        plotname_APdhw      = plotname_APdhw +'_smallCm.png'
        plotname_Nspikes_avgrms = plotname_Nspikes +'_smallCm_avgrms.png'
        plotname_APampl_avgrms  = plotname_APampl +'_smallCm_avgrms.png'
        plotname_APdhw_avgrms   = plotname_APdhw +'_smallCm_avgrms.png'
    else:
        plotname_Nspikes    = plotname_Nspikes +'.png'
        plotname_APampl     = plotname_APampl +'.png'
        plotname_APdhw      = plotname_APdhw +'.png'
        plotname_Nspikes_avgrms = plotname_Nspikes +'_avgrms.png'
        plotname_APampl_avgrms  = plotname_APampl +'_avgrms.png'
        plotname_APdhw_avgrms   = plotname_APdhw +'_avgrms.png'
    
    for i in range(Nmodels):
        testmodel   = testmodels[i]
        print('Model:',testmodel)
        Nspikes     = numpy.zeros(NCms)
        avg_AP_ampl = numpy.zeros(NCms)
        rms_AP_ampl = numpy.zeros(NCms)
        avg_AP_halfwidth = numpy.zeros(NCms)
        rms_AP_halfwidth = numpy.zeros(NCms)
        
        # Set names #_avgrms
        folder = 'Results/%i/IStim/current_idur'%testmodel+str(idur)+'_iamp'+str(iamp)+'/'
        outfilename_Nspikes = folder+'basps_idur%i_iamp'% (idur)+str(iamp) +'_Nspikes_vs_Cm'
        outfilename_APampl  = folder+'basps_idur%i_iamp'% (idur)+str(iamp) +'_APampl_vs_Cm'
        outfilename_APdhw   = folder+'basps_idur%i_iamp'% (idur)+str(iamp) +'_APdurhalfwidth_vs_Cm'
        if smallCms==True:
            outfilename_Nspikes = outfilename_Nspikes +'_smallCm.txt'
            outfilename_APampl  = outfilename_APampl +'_smallCm.txt'
            outfilename_APdhw   = outfilename_APdhw +'_smallCm.txt'
        else:
            outfilename_Nspikes = outfilename_Nspikes +'.txt'
            outfilename_APampl  = outfilename_APampl +'.txt'
            outfilename_APdhw   = outfilename_APdhw +'.txt'
        
        # make files
        outfile_Nspikes = open(outfilename_Nspikes,'w')
        outfile_APampl  = open(outfilename_APampl,'w')
        outfile_APdhw   = open(outfilename_APdhw,'w')
        for j in range(NCms):
            print('Step ', j+1, ' of', NCms)
            filename = folder+'basps_cm'+str(cms[j])+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_V.txt' 
            Nspikes[j], avg_AP_ampl[j], rms_AP_ampl[j], avg_AP_halfwidth[j], rms_AP_halfwidth[j], Nampl, Ndur = main(filename,idelay,idur)
            Nspikes_avg[j]           += Nspikes[j]
            Nampl_avg[j]             += Nampl
            Ndur_avg[j]              += Ndur
            AP_ampl_avg[j]           += avg_AP_ampl[j]*Nampl     # Weight it
            AP_halfwidth_avg[j]      += avg_AP_halfwidth[j]*Ndur # Weight it
            Nspikes_all[j,i]          = Nspikes[j]
            Nampl_all[j,i]            = Nampl
            Ndur_all[j,i]             = Ndur
            AP_ampl_rms[j]           += (rms_AP_ampl[j]*Nampl)**2
            AP_halfwidth_rms[j]      += (rms_AP_halfwidth[j]*Ndur)**2
            AP_ampl_all[j,i]          = avg_AP_ampl[j]
            AP_halfwidth_all[j,i]     = avg_AP_halfwidth[j]
            AP_ampl_rms_all[j,i]      = rms_AP_ampl[j]
            AP_halfwidth_rms_all[j,i] = rms_AP_halfwidth[j]
            outfile_Nspikes.write('%.5f %i\n' % (cms[j],Nspikes[j]))
            outfile_APampl.write('%.5f %.10f %.10f\n' % (cms[j],avg_AP_ampl[j],rms_AP_ampl[j]))
            outfile_APdhw.write('%.5f %.10f %.10f\n' % (cms[j],avg_AP_halfwidth[j],rms_AP_halfwidth[j]))
        Nspikes_bymodel.append(Nspikes)
        avg_AP_ampl_bymodel.append(avg_AP_ampl)
        rms_AP_ampl_bymodel.append(rms_AP_ampl)
        avg_AP_halfwidth_bymodel.append(avg_AP_halfwidth)
        rms_AP_halfwidth_bymodel.append(rms_AP_halfwidth)
        outfile_Nspikes.close()
        outfile_APampl.close()
        outfile_APdhw.close()


    for i in range(Ncms):
        AP_ampl_avg[i]       /= Nampl_avg[i] # Divide by number of peaks
        AP_halfwidth_avg[i]  /= Ndur_avg[i]
        AP_ampl_rms[i]       /= Nampl_avg[i]**2
        AP_halfwidth_rms[i]  /= Ndur_avg[i]**2
        AP_ampl_rms[i]        = np.sqrt(AP_ampl_rms[i])
        AP_halfwidth_rms[i]   = np.sqrt(AP_halfwidth_rms[i])
        Nspikes_avg[i] /= Nmodels       # Get average number of peaks
        Nampl_avg[i]  /= Nmodels       # Get average number of amplitudes ## Don't think I need these two
        Ndur_avg[i]   /= Nmodels       # Get average number of durations  ## Don't think I need these two
        Nspikes_temp = Nspikes_all[i,:]
        for j in range(Nmodels):
            Nspikes_rms[i] += (Nspikes_avg[i]-Nspikes_temp[j])**2
        Nspikes_rms[i] = np.sqrt(Nspikes_rms[i]/(Nmodels-1))
    
    ############################ Dendritepropagation #############################################
    idur = 2
    iamp = 1.0
    v_init = -65
    cms_dp = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
    Ncms_dp = len(cms_dp)
    propvels_avg = np.zeros(Ncms_dp)
    propvels_rms = np.zeros(Ncms_dp)
    propvels_all = np.zeros((Ncms_dp,Nmodels))
    i=0
    for testmodel in testmodels:
        Radp = 150
        if testmodel==488462965:
            Radp = 29
        elif testmodel==478513407:
            Radp = 100
        elif testmodel==478513437:
            Radp = 352
        dendfolder = 'Results/%i/IStim/current_idur'%testmodel+str(idur)+'_iamp'+str(iamp)+'/dendritepropagation/'
        filename_cms = dendfolder+'basps_varycm_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Radp+str(v_init)+'_propvel_vs_cm_smallCm.txt'
        
        propvels_temp = []
        
        file_cms = open(filename_cms,'r')
        lines = file_cms.readlines()
    
        j=0
        for line in lines:
            words = line.split()
            if len(words)>1:
                propvels_all[j,i] = float(words[1])
            j+=1
        i+=1
        file_cms.close()
    
    for i in range(Ncms_dp):
        propvels_avg[i],propvels_rms[i] = avg_and_rms(propvels_all[i,:])


    ########### Plot, rms and avg. ######################################################
    plottitletext = 'soma'
    print('cms:',cms)
    print('Nspikes_avg:',Nspikes_avg)
    print('len(cms):',len(cms))
    print('len(Nspikes_avg):',len(Nspikes_avg))
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(17,15))
    ax1.errorbar(cms, Nspikes_avg, yerr=Nspikes_rms, capsize=2)
    ax1.set(ylabel=r'Number of spikes')
    ax1.set_title(r'Number of spikes vs %s capacitance' % plottitletext)
    
    ax2.errorbar(cms, AP_ampl_avg, yerr=AP_ampl_rms, capsize=2)
    ax2.set(ylabel=r'Amplitude [mV]')
    ax2.set_title(r'Amplitude vs %s capacitance' % plottitletext)
    
    ax3.errorbar(cms, AP_halfwidth_avg, yerr=AP_halfwidth_rms, capsize=2)
    ax3.set(xlabel=r'$C_{m}$ of %s [$\mu$F/cm$^2$]' % plottitletext,ylabel=r'AP width at half amplitude [ms]')
    ax3.set_title(r'AP width at half amplitude vs %s capacitance' % plottitletext)
    
    ax4.errorbar(cms_dp, propvels_avg, yerr=propvels_rms, capsize=2)
    ax4.set(xlabel=r'$C_{m}$ of %s [$\mu$F/cm$^2$]' % plottitletext,ylabel=r'AP width at half amplitude [ms]')
    ax4.set_title(r'AP width at half amplitude vs %s capacitance' % plottitletext)
    
    plt.savefig(avgandrmsname)
    
    ###################### Plot, every model ############################################
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(17,15))
    for i in range(Nmodels):
        ax1.plot(cms,Nspikes_bymodel[i],'-o',label='%i' % testmodels[i])
    ax1.set(ylabel=r'$N_{spikes}$')
    ax1.set_title(r'Number of spikes vs soma capacitance')
    ax1.legend(loc='lower right')
        
    for i in range(Nmodels):
        ax2.errorbar(cms,avg_AP_ampl_bymodel[i], yerr=rms_AP_ampl_bymodel[i], capsize=2,label='%i' % testmodels[i])
    ax2.set(ylabel='Spike amplitude [mV]')
    ax2.set_title(r'AP amplitude vs soma capacitance')
    
    for i in range(Nmodels):
        ax3.errorbar(cms,avg_AP_halfwidth_bymodel[i], yerr=rms_AP_halfwidth_bymodel[i], capsize=2,label='%i' % testmodels[i])
    ax3.set(xlabel=r'$C_{m} $[$\mu$ F/cm$^2$]',ylabel=r'AP duration at half width [ms]')
    ax3.set_title(r'AP duration at half width vs soma capacitance')
    
    for i in range(Nmodels):
        ax4.plot(cms_dp,propvels_all[:,i], '-o', label='%i' % testmodels[i])
    ax4.set(xlabel=r'$C_{m} $[$\mu$ F/cm$^2$]',ylabel=r'Signal velocity [m/s]')
    ax4.set_title(r'Propagation velocity along dendrite vs soma capacitance')
    
    plt.savefig(everymodelname)

    