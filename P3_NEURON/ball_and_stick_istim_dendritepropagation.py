import numpy as np
import matplotlib.pyplot as plt
from neuron import h
import LFPy

def run_sim(cm=1.0,idur=1.0,iamp=1.0,idelay=1.0,v_init=-65.,tstart=0.,tstop=50.,Ra=150):

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
            'idx': stim_idx,
            'record_current': True,
            'pptype': 'IClamp',
            'amp': iamp,
            'dur': idur,
            'delay': idelay,
        }

    stimulus = LFPy.StimIntElectrode(cell, **stim_params)
    cell.simulate(rec_vmem=True, rec_imem=True)
    t, v = cell.tvec.copy(), cell.vmem[0].copy()
    Nt = len(t)
    
    idxs = cell.get_idx() # Maybe I should use this for picking sections
    #Print to file # Simplified because I will probably keep the sectioning, if not the entire morphology
    segmentlist = np.array([0,1,13,25,37,49]) # This is fine
    plotsegment = np.array([0,1,25,49]) # This is fine
    Nidx   = len(segmentlist)
    folder = 'Results/IStim/current_idur'+str(idur)+'_iamp'+str(iamp)+'/dendritepropagation/'
    plotname = folder + 'bas_cm'+str(cm)+'_idur%.1f_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_diffsegm.png'
    for i in range(Nidx):
        filename = folder+'bas_cm'+str(cm)+'_idur%.1f_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_V_seg%i.txt' %i # For easy access
        file = open(filename,'w')
        vsec = cell.vmem[segmentlist[i]].copy()
        for j in range(Nt):
            file.write('%.16f %.16f\n' % (t[j],vsec[j]))
        file.close()
    
    plt.figure(figsize=(6,5))
    # Plotting
    for i in range(len(plotsegment)):
       plt.plot(t, cell.vmem[plotsegment[i]].copy(), label='Segment %i' % plotsegment[i])
    plt.xlabel('Time [ms]')
    plt.ylabel('Voltage [mV]')
    plt.title(r'Voltage vs time for different cell segments (current input), $c_m$=%.1f' %cm)
    plt.legend()
    plt.savefig(plotname) 
             
    cell.__del__()
    return t, v # Do I ACTUALLY need to return something, though?

if __name__ == '__main__':
    cms = [0.5,1,2,3,4,4.5,5,5.5,6,7,8,9,10]
    iamp = 1.0
    idur = 1.0
    idelay  = 1.0
    afteri  = 10.0
    tstart  = 0.
    tstop_i = 50.#idur+afteri+idelay
    v_init  = -65
    Ra      = 150
    '''
    outfolder = 'Results/IStim/'
    currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    outfolder = outfolder+currentfolder
    plotname = outfolder+'bas_varycm_idur%.1f_iamp'%idur+str(iamp)+'_Ra%i_vinit'%Ra+str(v_init)+'_V.png'
    ''' 
    for cm in cms:
        t, v = run_sim(cm,idur,iamp,idelay,v_init,tstart,tstop_i,Ra)