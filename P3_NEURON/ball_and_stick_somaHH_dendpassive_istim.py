import numpy as np
import matplotlib.pyplot as plt
from neuron import h
import LFPy

def run_sim(cm=1.0,idur=1.0,iamp=1.0,idelay=1.0,v_init=-65.,tstart=0.,tstop=50.,Ra=150):

    dt = 2**-6 # 4
    cell_params = {          # various cell parameters,
                'morphology': "ballandstick.hoc",
                'v_init' : v_init,    # initial crossmembrane potential
                'cm': cm,
                'Ra': Ra,
                'passive': False,
                'passive_parameters': {'g_pas': 1./30000, 'e_pas': -65},
                'nsegs_method' : "lambda_f",
                'lambda_f': 500,
                'dt' : dt,   # [ms] dt's should be in powers of 2
                'tstart' : tstart,
                'tstop' : tstop,
            }


    cell = LFPy.Cell(**cell_params)
    h("soma insert hh")  # Inserting Hodgkin-Huxley mechanisms in soma
    h("dend insert pas") # Inserting passive mechanism in dendrite
    
    for sec in h.allsec():
        sectype = sec.name().split("[")[0]
        if sectype=='dend':
            sec.e_pas = -54.3
            sec.g_pas = 0.0003        

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
    cms = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0]
    iamp = 0.3
    idur = 1000
    idelay  = 100.0
    afteri  = 10.0
    tstart  = -100.
    tstop_i = idur+afteri+idelay
    v_init  = -65
    Ra      = 150
    outfolder = 'Results/IStim/'
    currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    outfolder = outfolder+currentfolder 
    for cm in cms:
        t, v = run_sim(cm,idur,iamp,idelay,v_init,tstart,tstop_i,Ra)
        outfilename = outfolder+'basHHdendpass_cm'+str(cm)+'_idur%.1f_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_V.txt' 
        outfile = open(outfilename,'w')
        # Write to 'Results/IStim/'

        Nt   = len(t)
        for i in range(Nt):
            outfile.write('%.16f %.16f\n' % (t[i],v[i]))
        outfile.close()    
        plotname = outfolder+'basHHdendpass_cm'+str(cm)+'_idur%.1f_iamp'%idur+str(iamp)+'_Ra%i_vinit'%Ra+str(v_init)+'_V.png'
        plt.figure(figsize=(6,5))
        plt.plot(t, v, label="c$_m$= {:1.2f} µF/cm²".format(cm))
        plt.xlabel('Time [ms]')
        plt.ylabel('Voltage [mV]')
        plt.title(r'Voltage vs time for different cell capacitances $c_m$ (current input)')
        plt.legend()
        plt.savefig(plotname)