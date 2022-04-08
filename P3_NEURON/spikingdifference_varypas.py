import numpy
import matplotlib.pyplot as plt

def reldiffs_one(fdefault,fother):
    rdout  = 0
    rdout2 = 0
    if fdefault!=0 and fother!=0:
        rdout = ((fdefault-fother)/float(fdefault))
        rdout2 = ((fdefault-fother)/float(fother))
    return rdout, rdout2


def avg_and_rms(x):
    N = len(x)
    avgx = numpy.mean(x)
    rmsx = 0
    for i in range(N):
        rmsx += (avgx-x[i])**2
    rmsx = numpy.sqrt(rmsx/(N-1))
    return avgx,rmsx

def manual(filename,idelay,idur,spikedurat):
    """Manual approach"""

    # Use numpy to read the trace data from the txt file
    data = numpy.loadtxt(filename)
    
    # Time is the first column
    time = data[:, 0]
    # Voltage is the second column
    voltage = data[:, 1]
    
    vmax = max(voltage) 
    vmin = min(voltage) 
    deltav = vmax-vmin
    vthr  = -20   # If there is a peak above this value, we count it # Old: vmax-0.15*deltav
    vprev = vthr-40 # A peak never kicks in at initiation, even if I change vthr
    durthr = spikedurat # Height at which we measure the duration. Need to be calibrated by peaks? # Was -20
    Npeaks = 0
    peaktimes = []
    timestart = idur/2.+idelay
    for i in range (1,len(voltage)-1):  
        if voltage[i-1]<voltage[i] and voltage[i+1]<voltage[i] and voltage[i]>vthr and time[i]>timestart:
            peaktimes.append(time[i])
            Npeaks+=1
        
        if time[i]>idur+idelay:
            break

    isi = []
    for i in range(1,len(peaktimes)):
        isi.append(peaktimes[i]-peaktimes[i-1])
    isi_avg, isi_rms = avg_and_rms(isi)
    
    ## Checking if we've got consistent firing.
    # Goals:
    # Checking if there's no firing in the last half of the stim. interval
    # Making sure that there is no flatlining: check for firing within 5*ISI (5 is chosen a bit at random here)
    if len(peaktimes)>0:
        if peaktimes[-1]<=(idur/2.+idelay): 
            Npeaks=0
            print('No consistent firing')
        if peaktimes[-1]<=(idur+idelay-5*isi_avg): # This should cover the first test, but better safe than sorry
            Npeaks=0
            print('No consistent firing')
        else:
            print('The neuron is firing!')
    
    return Npeaks

if __name__ == '__main__':
    testmodels = [478513437, 488462965, 478513407] 
    cm         = 1.0
    spikedurat = -40
    idur       = 2000 #100 # ms
    idelay     = 100
    v_init     = -83.7 # mV
    Ra         = 100
    somasize   = 10 # 15 # 
    dendlen    = 1000
    denddiam   = 1
    nsegments  = 200 
    
    varymech = 'leak' # 
    varyE_bool = True
    namestringfirst = ''
    varyE = [[0,27],[0,120],[0,13.5]]# Old: [[0,23.5,30],[0,103],[0,12,15]]
    iamps = [0.7,0.6,0.4]
    namestringfirst = namestringfirst + 'Epasplus'
    varyg = 'None'
    
    varylist = [] # Should be redundant
    plotstring = '_vary'
    varylist = varyE
    plotstring = plotstring + 'E'
      
    if varymech=='Na':
        folderstring = 'VaryNa/' 
        plotstring = plotstring + '_Na'
    elif varymech=='leak':
        folderstring = 'VaryLeak/'
        plotstring = plotstring + '_leak'
    elif varymech=='K':
        folderstring = 'VaryK/'
        plotstring = plotstring + '_K'
    
    ###########################################################
    
    NE = len(varyE)
    zeroind = 3 # Default
    for i in range(NE):
        if varyE[i]==0:
            zeroind = i
    
    Nspikes = numpy.zeros(NE)
    
    outfolder_all = 'figures/Comparemodels/Vrest/'
    plotname_all  = outfolder_all+'%s_Allen_cmfall'%namestringfirst+str(cm)+'_idur%i_'% idur+'_manual_Vrest_vs_Epas.png'
    outfilename_all  = outfolder_all+'%s_Allen_cmfall'%namestringfirst+str(cm)+'_idur%i_'% idur+'_manual_Vrestdiff_vs_Epas.txt'
    print('varyE:',varyE)
    # Set names
    Vrest_avg_all = []
    Vrest_rms_all = []
    k = 0
    for testmodel in testmodels:
        varyEthis = varyE[k]
        iamp = iamps[k]
        NE2 = len(varyEthis)
        if testmodel==480633479:
            v_init = -96.8#-83.7#-90#-86.5# # Have a test here too
        elif testmodel==496497595:
            v_init = -86.5
        elif testmodel==488462965:
            v_init = -86.5 # Maybe I should have changed this...
        elif testmodel==497230463:
            v_init = -90
        elif testmodel==497233075:
            v_init = -90
        elif testmodel==478513437:
            v_init = -86.8
        elif testmodel==478513407:
            v_init = -83.7
        elif testmodel==497233271:
            v_init = -90
        elif testmodel==489931686:
            v_init = -95.7
        elif testmodel==485694403:
            v_init = -88.8
        outfolder = 'figures/%i/' % testmodel
        outfilename = outfolder+'%s_%i_cmfall'%(namestringfirst,testmodel)+str(cm)+'_idur%i_'% idur+'_manual_Nspikes_vs_Epas.txt'
        outfilename_diff = outfolder+'%s_%i_cmfall'%(namestringfirst,testmodel)+str(cm)+'_idur%i_'% idur+'_manual_Nspikesdiff_vs_Epas.txt'
        plotname  = outfolder+'%s_%i_cmfall'%(namestringfirst,testmodel)+str(cm)+'_idur%i_'% idur+'_manual_Nspikes_vs_Epas.png'
        plotname_diff  = outfolder+'%s_%i_cmfall'%(namestringfirst,testmodel)+str(cm)+'_idur%i_'% idur+'_manual_Nspikesdiff_vs_Epas.png'
        # make files
        outfile = open(outfilename,'w')
        outfile_diff = open(outfilename_diff,'w')
        thr_reached=False
        Nspikes_these = numpy.zeros(NE)
        for j in range(NE2):
            E = varyEthis[j]
            namestring = namestringfirst+str(E)
            infolder = 'figures/%i/current_idur%i_iamp' % (testmodel,idur) + str(iamp)+'/'
            print('Step ', j+1, ' of', NE)
            filename = infolder+namestring+"_changecmf" + str(cm) + "_everywhere_vinit"+str(v_init)+"_addedRa.txt"
            Nspikes_these[j] = manual(filename,idelay,idur,spikedurat) 
            outfile.write('%.2f %i\n' % (E,Nspikes_these[j])) 
        outfile.close()
        # Rel diff.
        Nspikes0 = Nspikes_these[0]
        for j in range(1,NE2):
            diff1, diff2 = reldiffs_one(Nspikes0,Nspikes_these[j])
            outfile_diff.write('%.2f %.5e\n' % (varyEthis[j],diff1))
        outfile_diff.close()
        k+=1