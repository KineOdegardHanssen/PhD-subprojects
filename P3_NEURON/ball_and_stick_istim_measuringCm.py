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
    
    # Try to get information on compartments:
    # print out section information: # Works even though I do everything through LFPy
    if cm==1:
        for sec in h.allsec():
            h.psection()
        idxs = cell.get_idx(section="dend")
        print('idxs:',idxs)
    
    cell.__del__()
    return t, v


if __name__ == '__main__':
    cms = [0.5,1,2,3,4,4.5,5,5.5,6,7,8,9,10]
    iamp = -0.5
    idur = 1000.0
    idelay  = 100.0
    afteri  = 100.0
    tstart  = 0.
    tstop_i = idur+afteri+idelay
    v_init  = -65
    Ra      = 150
    outfolder = 'Results/IStim/'
    currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    outfolder = outfolder+currentfolder
    plotname = outfolder+'bas_varycm_idur%.1f_iamp'%idur+str(iamp)+'_Ra%i_vinit'%Ra+str(v_init)+'_V.png' 
    for cm in cms:
        print('On cm=',cm)
        t, v = run_sim(cm,idur,iamp,idelay,v_init,tstart,tstop_i,Ra)
        outfilename = outfolder+'bas_cm'+str(cm)+'_idur%.1f_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_V.txt' 
        outfile = open(outfilename,'w')
        # Write to 'Results/IStim/'

        Nt   = len(t)
        for i in range(Nt):
            outfile.write('%.16f %.16f\n' % (t[i],v[i]))
        outfile.close()    
        plt.plot(t, v, label="c$_m$= {:1.2f} µF/cm²".format(cm))
    plt.xlabel('Time [ms]')
    plt.ylabel('Voltage [mV]')
    plt.title(r'Voltage vs time for different cell capacitances $c_m$ (current input)')
    plt.legend()
    plt.savefig(plotname)