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
    print('data:',data)

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
    print('current2:',current2)
    # Set the 'I' (=current) key of the trace # Hope this helps # Or do I need the stimulus current?
    trace1['I'] = current2 #!!!!
    print('Time:',time)
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
    capacitance = capacitance*10**12*100 # in pF
    
    # Do I need this?
    # Stim Current sample index for detecting current of stimulation step
    #stim = trace1['I'][idelay+idur]
    
    return capacitance

if __name__ == '__main__':
    testmodel = 496497595
    idur = 100 # ms
    iamp = -0.5 # nA
    Ra   = 150
    idelay = 100
    v_init = -86.5
    d = 20    # Soma diameter and length
    A = numpy.pi*d**2 # Cell area
    
    cms = [0.5,1,2,3,4,4.5,5,5.5,6,7,8,9,10]
    Ncm = len(cms)
    capacitances = numpy.zeros(Ncm)
    
    folder = 'Results/%i/IStim/current_idur'%testmodel+str(idur)+'_iamp'+str(iamp)+'/'
    outfilename_all = folder +'somaonly_%i_varycm_idur%i_iamp'%(testmodel,idur)+str(iamp)+'_Ra%i_vinit'%Ra+str(v_init)+'_Cm.txt' 
    plotname =  folder +'somaonly_%i_varycm_idur%i_iamp'%(testmodel,idur)+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_Cm.png' 
    outname_all = open(outfilename_all,'w')
    for i in range(Ncm):
        cm = cms[i]
        filename = folder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_V.txt' 
        outfilename = folder +'somaonly_%i_cm'%testmodel+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_Cm.txt' 
        capacitance = main(filename,idelay,idur,iamp)/A ### Per area.
        capacitances[i] = capacitance
        outfile = open(outfilename,'w')
        outfile.write('%.12f' % capacitance)
        outname_all.write('%.1f %.12f\n' % (cm,capacitance))
        outfile.close()
        print('total capacitance:', capacitance)
outname_all.close()

plt.figure(figsize=(6,5))
plt.plot(cms, capacitances,'-o')
plt.xlabel(r'Parameter $C_{m}$ [$\mu$ F/cm$^2$]')
plt.ylabel(r'Output $C_{m}$ [$\mu$ F/cm$^2$]')
plt.title(r'Input $C_{m}$ vs measured $C_{m}$')
plt.savefig(plotname)