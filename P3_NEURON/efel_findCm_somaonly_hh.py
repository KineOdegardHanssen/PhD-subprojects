import efel
import numpy
import math
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
    
    Vmax = max(voltage)
    Vmin = min(voltage)
    tmax = max(time)
    tmin = min(time)
    V0 = voltage[0]
    dV = V0-Vmin
    dV1e = dV*math.exp(-1)
    V1e  = Vmin+dV1e
    
    # Will perform linear interpolation near the crossing with V1e.
    tbef = 0
    taft = 0
    Vbef = 0
    Vaft = 0
    for i in range(len(voltage)):
        if voltage[i]-V1e<0:
            tbef = time[i-1]
            taft = time[i]
            Vbef = voltage[i-1]
            Vaft = voltage[i]
            print('V:',voltage[i-1],'i-1:',i-1,'; t:', time[i-1])
            print('V:',voltage[i],'i:',i,'; t:', time[i])
            break
    
    a = (Vaft-Vbef)/(taft-tbef)
    t_interpolated = (V1e-Vbef+a*tbef)/a-idelay
    
    print('tbef',tbef)
    print('taft',taft)
    print('t_interpolated',t_interpolated)
    
    tau3 = t_interpolated
    tau4 = time_constant
    
    print('dV+Vmax:',dV+Vmax)
    testit = True
    if testit==True:
        N = 100*idur
        timecut = numpy.linspace(0,idur,N)
        fit3 = numpy.zeros(N)
        fit4 = numpy.zeros(N)
        for i in range(N):
            fit3[i] = dV*math.exp(-timecut[i]/tau3)+Vmin
            fit4[i] = dV*math.exp(-timecut[i]/tau4)+Vmin
        timecut+=idelay

        plottauy = numpy.array([Vmin-10,V0+10])
        plottaux = numpy.array([time_constant+idelay,time_constant+idelay])
        plotVex  = numpy.array([tmin,tmax])
        plotVey  = numpy.array([V1e,V1e])
        
        plt.figure(figsize=(6,5))
        plt.plot(time,voltage)
        plt.plot(timecut,fit3,'--',label='Manual fit')
        plt.plot(timecut,fit4,':',label='eFEL fit')
        #plt.plot(plottaux,plottauy,'c--')
        plt.plot(plotVex,plotVey,'c--')
        plt.xlabel('t')
        plt.ylabel('V')
        plt.axis([idelay-1,idelay+6*time_constant[0],Vmin-10,V0+10])
        plt.legend(loc='upper right')
        plt.show()
    
    print('time constant: eFEL:',time_constant[0], '; manual:', t_interpolated)
    print('time constant, eFEL-manual:', time_constant-t_interpolated)
    return time_constant[0], capacitance[0]

if __name__ == '__main__':
    idur = 1000 # ms
    iamp = -0.1 # nA
    Ra   = 100 #150
    idelay = 10 #100
    v_init = -70 #-86.5
    d = 10    # Soma diameter and length
    A = numpy.pi*d**2 # Cell area
    
    cms = [0.5,0.8,1.0,1.2,1.5]
    Ncm = len(cms)
    capacitances = numpy.zeros(Ncm)
    time_constants = numpy.zeros(Ncm)
    
    folder = 'Results/IStim/Soma%i/current_idur'%d+str(idur)+'_iamp'+str(iamp)+'/'
    outfilename_all = folder +'somaonly_varycm_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit'%Ra+str(v_init)+'_Cm.txt' 
    outfilename_all_tau = folder +'somaonly_varycm_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit'%Ra+str(v_init)+'_tau.txt' 
    plotname =  folder +'somaonly_varycm_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_Cm.png' 
    outname_all = open(outfilename_all,'w')
    outname_all_tau = open(outfilename_all_tau,'w')
    for i in range(Ncm):
        cm = cms[i]
        print('--------------------------------')
        print('Cm:',cm)
        filename = folder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_V.txt' 
        outfilename = folder +'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_Cm.txt' 
        time_constant, capacitance = main(filename,idelay,idur,iamp) ### Per area.
        capacitance = capacitance/A
        capacitances[i] = capacitance
        time_constants[i] = time_constant
        outfile = open(outfilename,'w')
        outfile.write('%.12f' % capacitance)
        outname_all.write('%.1f %.12f\n' % (cm,capacitance))
        outname_all_tau.write('%.1f %.12f\n' % (cm,time_constant))
        outfile.close()
        print('total capacitance:', capacitance)
outname_all.close()

plt.figure(figsize=(6,5))
plt.plot(cms, capacitances,'-o')
plt.xlabel(r'Parameter $C_{m}$ [$\mu$ F/cm$^2$]')
plt.ylabel(r'Output $C_{m}$ [$\mu$ F/cm$^2$]')
plt.title(r'Input $C_{m}$ vs measured $C_{m}$')
plt.savefig(plotname)

plt.figure(figsize=(6,5))
plt.plot(cms, capacitances,'-o')
plt.xlabel(r'Parameter $C_{m}$ [$\mu$ F/cm$^2$]')
plt.ylabel(r'Output $C_{m}$ [$\mu$ F/cm$^2$]')
plt.title(r'Input $C_{m}$ vs measured $C_{m}$')
plt.savefig(plotname)

plt.show()