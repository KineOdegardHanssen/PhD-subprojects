
import numpy as np
import matplotlib.pyplot as plt
from neuron import h
import LFPy

def run_sim(cm=1.0):

    dt = 2**-4
    cell_params = {          # various cell parameters,
                'morphology': "ballandstick.hoc",
                'v_init' : -65.,    # initial crossmembrane potential
                'cm': cm,
                'Ra': 150,
                'passive': True,
                'passive_parameters': {'g_pas': 1./30000, 'e_pas': -65},
                'nsegs_method' : "lambda_f",
                'lambda_f': 500,
                'dt' : dt,   # [ms] dt's should be in powers of 2
                'tstart' : 0.,
                'tstop' : 50.,
            }


    cell = LFPy.Cell(**cell_params)
    h("soma insert hh") # Inserting spiking mechanism

    stim_idx = 0
    stim_params = {
                 'idx' : stim_idx,
                 'record_current' : True,
                 'syntype' : 'Exp2Syn',
                 'tau1': 0.1,
                 'tau2': 0.2,
                 'weight' : 0.01,
                }

    syn = LFPy.Synapse(cell, **stim_params)
    syn.set_spike_times(np.array([10]))
    cell.simulate(rec_vmem=True, rec_imem=True)
    t, v = cell.tvec.copy(), cell.vmem[0].copy()
    cell.__del__()
    return t, v


if __name__ == '__main__':

    cms = [0.5, 1, 2, 10]
    for cm in cms:
        t, v = run_sim(cm)

        plt.plot(t, v, label="c$_m$= {:1.2f} µF/cm²".format(cm))
    plt.legend()
    plt.show()