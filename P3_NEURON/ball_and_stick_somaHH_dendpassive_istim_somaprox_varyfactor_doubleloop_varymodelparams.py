import numpy as np
import matplotlib.pyplot as plt
from neuron import h
import LFPy

def run_sim(varymech, varyE, varyg, somasize, cm=1.0, cmfac=1.0, idur=1.0, iamp=1.0, idelay=1.0, v_init=-65., tstart=-100., tstop=50., Ra=150, nsegments=200, dtexp=-8):

    dt = 2**dtexp
    cell_params = {          # various cell parameters,
                'morphology': "ballandstick.hoc",
                'v_init' : v_init,    # initial crossmembrane potential
                'cm': cm,
                'Ra': Ra,
                'passive': False,
                'passive_parameters': {'g_pas': 1./30000, 'e_pas': -65},
                'nsegs_method' : "lambda_f",
                'lambda_f': 200,
                'dt' : dt,   # [ms] dt's should be in powers of 2
                'tstart' : tstart,
                'tstop' : tstop,
            }


    cell = LFPy.Cell(**cell_params)
    h("soma insert hh")  # Inserting Hodgkin-Huxley mechanisms in soma
    h("dend insert pas") # Inserting passive mechanism in dendrite
    
    # Changing mechanisms
    # AND
    # Changing cm in proximal regions
    cmthis = cm*cmfac
    pnncutoff = 3.5*somasize
    for sec in h.allsec():
        sec.Ra = Ra
        sectype = sec.name().split("[")[0]
        
        if sectype=="soma": # Works
            # First: Change mechanisms
            if varymech=='Na':
                if varyE!='None':
                    sec.ena = varyE
                if varyg!='None':
                    sec.gnabar_hh = varyg
            elif varymech=='pas':
                if varyE!='None':
                    sec.el_hh = varyE
                if varyg!='None':
                    sec.gl_hh = varyg
            elif varymech=='K':
                if varyE!='None':
                    sec.ek = varyE
                if varyg!='None':
                    sec.gkbar_hh = varyg
            exec("sec.cm = {}".format(cmthis))
            sec.L = somasize
            sec.diam = somasize
            somasec  = sec
        if sectype=="dend":
            # First: Change passive mechanisms
            sec.e_pas = -54.3 
            sec.g_pas = 0.0003     
            # Then: Change geometry
            sec.L = dendlen
            sec.diam = denddiam
            sec.nseg = nsegments
            section = sec
            dist = h.distance(somasec(0.5),section(1)) # Measuring distance to end of section
            if dist<=pnncutoff:
                exec("sec.cm = {}".format(cmthis))
            else:
                for seg in sec:
                    if h.distance(somasec(0.5),seg)<=pnncutoff:
                        exec("seg.cm = {}".format(cmthis))
    
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
    
    cell.__del__()
    return t, v


if __name__ == '__main__':
    varymech = 'Na' # 'K' # 'leak'
    varyE_bool = True
    varyE = [40,50,60,70]
    varyg = 'None'
    
    varylist = [] # Should be redundant
    plotstring = '_vary'
    if varyE_bool==True:
        varylist = varyE
        plotstring = plotstring + 'E'
    else:
        varylist = varyg
        plotstring = plotstring + 'g'
    Nvary    = len(varylist)
      
    if varymech=='Na':
        folderstring = 'VaryNa/' 
        plotstring = plotstring + '_Na'
    elif varymech=='leak':
        folderstring = 'VaryLeak/'
        plotstring = plotstring + '_leak'
    elif varymech=='K':
        folderstring = 'VaryK/'
        plotstring = plotstring + '_K'
    
    dtexp     = -8
    cm        = 1.0 
    cmfac     = 1.0
    iamps     = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5] # Maybe I should divide this in half... # Wanted to do linspace, but ran into troubles before...
    idur      = 1000
    idelay    = 100.0
    afteri    = 10.0
    tstart    = -100.
    tstop_i   = idur+afteri+idelay
    v_init    = -65
    Ra        = 100
    somasize  = 10
    dendlen   = 1000
    denddiam  = 1
    nsegments = 200
    outfolder_base = 'Results/IStim/Soma%i/dendlen%i/denddiam%i/' % (somasize,dendlen,denddiam)+ folderstring
    for iamp in iamps:
        currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
        outfolder = outfolder_base+currentfolder
        plotname = outfolder+'basHHdpas_cmf'+str(cmfac)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit'%Ra+str(v_init)+plotstring+'_V.png'
        #plt.figure(figsize=(6,5))
        for elem in varylist:
            print('elem:',elem, '; current:', iamp)
            changestring =''
            if varyE_bool==True:
                varyE = elem
                changestring = changestring+'_E'+str(varyE)+'_gdflt'
            else:
                varyg = elem
                changestring = changestring+'_Edefault_g'+str(varyg)
            t, v = run_sim(varymech, varyE, varyg, somasize, cm, cmfac, idur, iamp, idelay, v_init, tstart, tstop_i, Ra, nsegments, dtexp)
            outfilename = outfolder+'basHHdpas_cmf'+str(cmfac)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+changestring+'_V.txt' 
            outfile = open(outfilename,'w')
            # Write to 'Results/IStim/'
    
            Nt   = len(t)
            for i in range(Nt):
                outfile.write('%.16f %.16f\n' % (t[i],v[i]))
            outfile.close()    
            #plt.plot(t, v, label="c$_m$= {:1.2f} µF/cm²".format(cm))
        #plt.xlabel('Time [ms]')
        #plt.ylabel('Voltage [mV]')
        #plt.title(r'Voltage vs time for different cell capacitances $c_m$ (current input)')
        #plt.legend()
        #plt.savefig(plotname)