import os
from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import json
import neuron
import LFPy

def run_sim(varymech, varyE, varyg,somasize,dendlen,denddiam,nsegments,cm=1.0, idur=1.0, iamp=1.0, idelay=1.0, v_init=-86.5, tstart=0., tstop=50., testmodel=496497595, dtexp=-8):
    
    ### Copy-pasted, more or less, from the test_allan_cell_models-script
    testmodelname = 'neur_%i' % testmodel
    model_folder = 'cell_models/%s/' % testmodelname
    params = json.load(open(join(model_folder, "fit_parameters.json"), 'r'))
    
    e_pas = params["passive"][0]["e_pas"]
    celsius = params["conditions"][0]["celsius"]
    reversal_potentials = params["conditions"][0]["erev"]
    #v_init = params["conditions"][0]["v_init"] # This might be wrong, will set it by hand above
    active_mechs = params["genome"]
    Ra = params["passive"][0]["ra"]  # This was not here before
    neuron.h.celsius = celsius
    
    dt = 2.**dtexp #2.**-6 #-4
    cell_params = {          # various cell parameters,
                'morphology': "ballandstick.hoc",
                'v_init' : v_init,    # initial crossmembrane potential
                'cm': 1.0,
                'Ra': Ra,
                'passive': False,
                'nsegs_method' : "lambda_f",
                'lambda_f': 200.,
                'dt' : dt,   # [ms] dt's should be in powers of 2
                'tstart' : tstart,
                'tstop' : tstop,
            }


    cell = LFPy.Cell(**cell_params)
    
    pnncutoff = 3.5*somasize
    for sec in neuron.h.allsec():
        sec.insert("pas")
        sec.e_pas = e_pas
        sec.Ra = Ra
        sectype = sec.name().split("[")[0]
        
        for sec_dict in active_mechs:
            if sec_dict["section"] == sectype:
                if sectype=="soma": # Works
                    exec("sec.cm = {}".format(cm))
                    sec.L = somasize
                    sec.diam = somasize
                    somasec  = sec
                if sectype=="dend":
                    sec.L = dendlen
                    sec.diam = denddiam
                    sec.nseg = nsegments
                    section = sec
                    dist = neuron.h.distance(somasec(0.5),section(1)) # Measuring distance to end of section
                    if dist<=pnncutoff:
                        exec("sec.cm = {}".format(cm))
                    else:
                        for seg in sec:
                            if neuron.h.distance(somasec(0.5),seg)<=pnncutoff:
                                exec("seg.cm = {}".format(cm)) # 
                if not sec_dict["mechanism"] == "":
                    sec.insert(sec_dict["mechanism"])
                exec("sec.{} = {}".format(sec_dict["name"], sec_dict["value"]))

        for sec_dict in reversal_potentials:
            if sec_dict["section"] == sectype:
                # print(sectype, sec_dict)
                for key in sec_dict.keys():
                    if not key == "section":
                        exec("sec.{} = {}".format(key, sec_dict[key]))

    stim_idx = 0
    stim_params = {
            'idx': stim_idx,
            'record_current': True,
            'pptype': 'IClamp',
            'amp': iamp,
            'dur': idur,
            'delay': idelay,
        }
    
    # Vary properties:
    for sec in neuron.h.soma:
        if varymech=='NaV':
            if varyE!='None':
                sec.ena = varyE
            if varyg!='None':
                sec.gbar_NaV = varyg
        elif varymech=='pas':
            if varyE!='None':
                sec.e_pas = varyE
            if varyg!='None':
                sec.g_pas = varyg
        elif varymech=='Kd':
            if varyE!='None':
                sec.ek = varyE
            if varyg!='None':
                sec.gbar_Kd = varyg
        elif varymech=='Kv2like':
            if varyE!='None':
                sec.ek = varyE
            if varyg!='None':
                sec.gbar_Kv2like = varyg
        elif varymech=='Kv3_1':
            if varyE!='None':
                sec.ek = varyE
            if varyg!='None':
                sec.gbar_Kv3_1 = varyg
        elif varymech=='SK':
            if varyE!='None':
                sec.ek = varyE
            if varyg!='None':
                sec.gbar_SK = varyg
        elif varymech=='K_T':
            if varyE!='None':
                sec.ek = varyE
            if varyg!='None':
                sec.gbar_K_T = varyg
        elif varymech=='Im_v2':
            if varyE!='None':
                sec.ek = varyE
            if varyg!='None':
                sec.gbar_Im_v2 = varyg

    stimulus = LFPy.StimIntElectrode(cell, **stim_params)
    cell.simulate(rec_vmem=True, rec_imem=True)
    t, v = cell.tvec.copy(), cell.vmem[0].copy()
    
    # Try to get information on compartments:
    # print out section information: # Works even though I do everything through LFPy
    '''
    if cm==0.5:
        for sec in neuron.h.allsec():
            neuron.h.psection()
            for seg in sec:
                print('seg.cm:',seg.cm, '; seg.x:', seg.x, 'seg.diam:', seg.diam, 'seg.g_pas:',seg.g_pas, 'seg.e_pas:',seg.e_pas)
 
        #idxs = cell.get_idx(section="dend")
        #print('idxs:',idxs)
    '''
    
    cell.__del__()
    return t, v


if __name__ == '__main__':
    varymech = 'NaV' #'SK'#'Im_v2'#'Kv2like'#'Kd'#'pas' #
    varyE_bool = True
    varyE = [40,50,53,65,70]#[-60,-53,-30,-10,10,30,60] # Default...407Kd-107 Default478513407pas: -83.6528;  Default478513407Na: 53
    varyg = 'None' #0.09 # Default...407Kd2.94396e-010 # Default478513407pas: 0.000362109; Default478513407Na: 0.0409177
    
    varylist = [] # Should be redundant
    plotstring = '_vary'
    if varyE_bool==True:
        varylist = varyE
        plotstring = plotstring + 'E'
    else:
        varylist = varyg
        plotstring = plotstring + 'g'
    Nvary    = len(varylist)
      
    if varymech=='NaV':
        folderstring = 'VaryNa/' 
        plotstring = plotstring + '_NaV'
    elif varymech=='pas':
        folderstring = 'VaryPas/'
        plotstring = plotstring + '_Pas'
    elif varymech=='Kd':
        folderstring = 'VaryKd/'
        plotstring = plotstring + '_Kd'
    elif varymech=='Kv2like':
        folderstring = 'VaryKv2like/'
        plotstring = plotstring + '_Kv2like'
    elif varymech=='Kv3_1':
        folderstring = 'VaryKv3_1/'
        plotstring = plotstring + '_Kv3_1'
    elif varymech=='SK':
        folderstring = 'VarySK/'
        plotstring = plotstring + '_SK'
    elif varymech=='K_T':
        folderstring = 'VaryK_T/'
        plotstring = plotstring + '_K_T'
    elif varymech=='Im_v2':
        folderstring = 'VaryIm_v2/'
        plotstring = plotstring + '_Im_v2'
    dtexp = -8
    testmodel = 489931686 # 485694403 # 488462965 # 478513437 # 478513407 # 489931686 #   ### Weirdo:480633479#  ###AA:496497595
    smallCm = True
    cm = 1
    iamps = [0]#,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44] # Maybe I should divide this in half... # Wanted to do linspace, but ran into troubles before...
    idur = 1000
    idelay   = 100.0
    afteri   = 100.0
    tstart   = -500.
    tstop    = idur+afteri+idelay
    v_init   = -86.5
    somasize = 10 # 15 # 
    dendlen  = 1000
    denddiam = 2
    nsegments = 200 #dendlen # Not entirely sure about this variable. 1 mu m is probably unneccessarily small...
    outfolder_base = 'Results/%i/IStim/Soma%i/dendlen%i/denddiam'% (testmodel,somasize,dendlen)+str(denddiam)+'/'+ folderstring
    for iamp in iamps:
        print('iamp:',iamp)
        currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
        outfolder = outfolder_base+currentfolder
        plotname = outfolder+'baspv_varycm_idur%i_iamp'%idur+str(iamp)+'_vinit'+str(v_init)+plotstring+'_sprx_V.png'
        plt.figure(figsize=(6,5))
        for elem in varylist:
            print('elem:',elem)
            changestring =''
            if varyE_bool==True:
                varyE = elem
                changestring =changestring+'_E'+str(varyE)+'_gdflt'
            else:
                varyg = elem
                changestring =changestring+'_Edefault_g'+str(varyg)
            t, v = run_sim(varymech, varyE, varyg, somasize,dendlen,denddiam,nsegments,cm, idur, iamp, idelay, v_init, tstart, tstop, testmodel,dtexp)
            outfilename = outfolder+'baspv_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_vinit'+str(v_init)+changestring+'_sprx_V.txt' 
            outfile = open(outfilename,'w')
            # Write to 'Results/IStim/'
    
            Nt   = len(t)
            for i in range(Nt):
                outfile.write('%.16f %.16f\n' % (t[i],v[i]))
            outfile.close()    
            plt.plot(t, v, label='param=%.2f' % elem)
        plt.xlabel('Time [ms]')
        plt.ylabel('Voltage [mV]')
        plt.title(r'Voltage vs time for different parameters, %s' % plotstring)
        plt.legend()
        plt.savefig(plotname)