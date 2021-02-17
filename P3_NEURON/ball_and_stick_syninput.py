import numpy as np
import matplotlib.pyplot as plt
from neuron import h
import LFPy

def run_sim(cm=1.0,tau1=0.1,tau2=0.2,weight=0.01,v_init=-65.,tstart=0.,tstop=50.,Ra=150):

    dt = 2**-4
    cell_params = {          # various cell parameters,
                'morphology': "ballandstick.hoc",
                'v_init' : v_init,    # initial crossmembrane potential
                'cm': cm,
                'Ra': Ra,
                'passive': True,
                'passive_parameters': {'g_pas': 1./30000, 'e_pas': -65},
                'nsegs_method' : "lambda_f",
                'lambda_f': 500,
                'dt' : dt,   # [ms] dt's should be in powers of 2
                'tstart' : tstart,
                'tstop' : tstop,
            }


    cell = LFPy.Cell(**cell_params)
    h("soma insert hh") # Inserting spiking mechanism

    stim_idx = 0
    stim_params = {
                 'idx' : stim_idx,
                 'record_current' : True,
                 'syntype' : 'Exp2Syn',
                 'tau1': tau1,
                 'tau2': tau2,
                 'weight' : weight,
                }

    syn = LFPy.Synapse(cell, **stim_params)
    syn.set_spike_times(np.array([10]))
    cell.simulate(rec_vmem=True, rec_imem=True)
    t, v = cell.tvec.copy(), cell.vmem[0].copy()
    cell.__del__()
    return t, v


if __name__ == '__main__':
    outfolder = 'Results/SynInput/'
    cms = [0.5, 1, 2, 10]
    tau1 = 0.1
    tau2 = 0.2
    weight = 0.01
    v_init =-65
    tstart = 0.
    tstop  = 50.
    Ra     = 150
    plotname = outfolder+'bas_varycm_tau1%.1f_tau2%.1f_weight%.2f_Ra%i_vinit'% (tau1,tau2,weight,Ra)+str(v_init)+'_V.png' 
    for cm in cms:
        t, v = run_sim(cm,tau1,tau2,weight,v_init,tstart,tstop,Ra)
        outfilename = outfolder+'bas_cm'+str(cm)+'_tau1%.1f_tau2%.1f_weight%.2f_Ra%i_vinit'% (tau1,tau2,weight,Ra)+str(v_init)+'_V.txt' 
        outfile = open(outfilename,'w')
        # Write to 'Results/SynInput/'

        Nt   = len(t)
        for i in range(Nt):
            outfile.write('%.16f %.16f\n' % (t[i],v[i]))
        outfile.close()    
        plt.plot(t, v, label="c$_m$= {:1.2f} µF/cm²".format(cm))
    plt.xlabel('Time [ms]')
    plt.ylabel('Voltage [mV]')
    plt.title(r'Voltage vs time for different cell capacitances $c_m$ (synaptic input)')
    plt.legend()
    plt.savefig(plotname)