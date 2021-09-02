import numpy as np
import matplotlib.pyplot as plt
from neuron import h
import LFPy

def run_sim(somasize,dendlen,denddiam,nsegments,cm=1.0,idur=1.0,iamp=1.0,idelay=1.0,v_init=-65.,tstart=0.,tstop=50.,Ra=100, dtexp=-8,gpas=0.0003,vpas=-65):

    dt = 2**dtexp #-8 #6 # 4
    cell_params = {          # various cell parameters,
                'morphology': "ballandstick.hoc",
                'v_init' : v_init,    # initial crossmembrane potential
                'cm': cm,
                'Ra': Ra,
                'passive': False,
                'passive_parameters': {'g_pas': 0.001, 'e_pas': -65},  # Rewrite this later on.
                'nsegs_method' : "lambda_f",
                'lambda_f': 500,
                'dt' : dt,   # [ms] dt's should be in powers of 2
                'tstart' : tstart,
                'tstop' : tstop,
            }


    cell = LFPy.Cell(**cell_params)
    
    # Set length
    for sec in h.allsec():
        sec.insert("pas")
        sectype = sec.name().split("[")[0]
        sec.cm = cm
        if sectype=="soma": # Works
            sec.L = somasize
            sec.diam = somasize
        if sectype=="dend":
            sec.L = dendlen
            sec.diam = denddiam
            sec.nseg = nsegments
        for seg in sec:
            seg.g_pas = gpas #0.0003: Standard value (HH)
            seg.e_pas = vpas

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
        #idxs = cell.get_idx(section="dend")
        #print('idxs:',idxs)
    
    cell.__del__()
    return t, v


if __name__ == '__main__':
    dtexp = -8 # -10
    cms = [0.1,0.5,0.75,0.88,0.94,1.0,1.25,1.5,2.0,3.0] 
    # cms = [0.01,0.1,0.5,0.8,1.0,1.2,1.5,2.0,3.0,5.0,7.0,10.0,15.0]
    iamp = -0.1 # 0.1 # 
    idur = 100
    idelay   = 10.0
    afteri   = 10.0
    tstart   = -100.
    tstop_i  = idur+afteri+idelay
    v_init   = -65
    Ra       = 100
    somasize = 10 # 15 # 
    dendlen  = 1
    denddiam = 20
    gpas     = 0.0003
    vpas     = -65
    nsegments = dendlen # Not entirely sure about this variable
    outfolder = 'Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/' 
    currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/dtexp%i/' % dtexp
    outfolder = outfolder+currentfolder 
    for cm in cms:
        t, v = run_sim(somasize,dendlen,denddiam,nsegments,cm,idur,iamp,idelay,v_init,tstart,tstop_i,Ra,dtexp,gpas,vpas)
        outfilename = outfolder+'baspass_cm'+str(cm)+'_idur%.1f_iamp'%idur+str(iamp)+'_Ra%i_gpas'%Ra+str(gpas)+'_vpas' +str(vpas)+'_V.txt' 
        outfile = open(outfilename,'w')
        # Write to 'Results/IStim/'

        Nt   = len(t)
        for i in range(Nt):
            outfile.write('%.16f %.16f\n' % (t[i],v[i]))
        outfile.close()    
        plotname = outfolder+'baspass_cm'+str(cm)+'_idur%.1f_iamp'%idur+str(iamp)+'_Ra%i_gpas'%Ra+str(gpas)+'_vpas' +str(vpas)+'_V.png' 
        plt.figure(figsize=(6,5))
        plt.plot(t, v, label="c$_m$= {:1.2f} µF/cm²".format(cm))
        plt.xlabel('Time [ms]')
        plt.ylabel('Voltage [mV]')
        plt.title(r'Voltage vs time for different cell capacitances $c_m$ (current input)')
        plt.legend()
        plt.savefig(plotname)
        
        if cm==1:
            for sec in h.allsec():
                h.psection()