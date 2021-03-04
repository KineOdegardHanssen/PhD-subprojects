import efel
import numpy
import matplotlib.pyplot as plt

## Mainly Sverre's stuff.
## Needed to modify to fit my simulations, though
## Setting I is a bit hacky.

def main(filename,idelay,idur,iamp):
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
    
    Nt = len(time)
    current2 = []
    current = numpy.zeros(Nt) # Or set this through LFPy?
    for i in range(Nt):
        if time[i]>idelay and time[i]<(idelay+idur):
             current[i] = iamp
             current2.append(iamp)
    current2 = numpy.array(current2)
    
    # Set the 'I' (=current) key of the trace # Hope this helps # Or do I need the stimulus current?
    trace1['I'] = current2 #!!!!
    print('current:',current)
    print('numpy.sum(current):',numpy.sum(current))
    
    # Fix this:
    # Stim Current sample index for detecting current of stimulation step
    trace1['stimulus_current'] = current2 #iamp #current

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

    v_deflect = (efel.getFeatureValues(traces,['voltage_deflection']))
    trace1['voltage_deflection']= v_deflect[0]['voltage_deflection']
    
    input_res = (efel.getFeatureValues(traces,['ohmic_input_resistance']))
    print('input_res:',input_res)
    input_res = (input_res[0]['ohmic_input_resistance'])#*1000

    # Time constant, see Efel docs.
    time_constant = efel.getFeatureValues(traces,['time_constant'])
    time_constant = time_constant[0]['time_constant']
    print('time_constant:',time_constant)
    
    # Capacitance in Farads.
    capacitance = (time_constant*10**-3)/(input_res*10**6)
    capacitance = capacitance*10**12 # in pF
    capacitance = capacitance*100 # To get the answer in uF/cm2 when I divide by um2
    
    # Do I need this?
    # Stim Current sample index for detecting current of stimulation step
    #stim = trace1['I'][idelay+idur]
    
    return capacitance

if __name__ == '__main__':
    testmodel = 496497595 #488462965 #
    idur = 1000 # ms
    iamp = -0.5 # nA
    idelay = 100
    idur   = 1000
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
    
    # Changing values of membrane capacitance:
    #cm_soma = 3.0
    cm_dend = 2.0 #cm_soma # 15.0 #10.0
    #cm_axon = cm_soma #10.0
    
    folder = 'Allen_test_changecapacitance/figures/%i/current_idur%i_iamp' % (testmodel,idur)+str(iamp)+'/'
    filename = folder+'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) + '_vinit'+str(v_init)+'_addedRa.txt'
    outfilename = folder +'idur%i_iamp' % idur+str(iamp)+'_cmsoma' + str(cm_soma) + '_cmdend' + str(cm_dend) + '_cmaxon'+ str(cm_axon) + '_vinit'+str(v_init)+'_addedRa_Cmtotal.txt'
    capacitance = main(filename,idelay,idur,iamp)
    outfile = open(outfilename,'w')
    outfile.write('%.12f' % capacitance)
    outfile.close()
    print('total capacitance:', capacitance)