import numpy
import matplotlib.pyplot as plt

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
    
    # Find resting potential:
    N = len(time)
    Nstart = int(0.3*N)
    tVrest = time[Nstart:N]
    Vrest  = 0 # Or store in array? Don't want overflow.
    Vrest  = []
    Vrest_rms = 0 # Should stay 0 to, or stay very close to 0. More of a sanity check.
    for i in range(Nstart,N):
        Vrest.append(voltage[i]) # Could append time too...
    Nrest = len(Vrest)
    Vrest = numpy.array(Vrest)
    Vrest_avg = numpy.mean(Vrest)
    for i in range(Nrest):
        Vrest_rms+=(Vrest[i]-Vrest_avg)**2
    Vrest_rms=numpy.sqrt(Vrest_rms/(Nrest-1))
    
    '''
    plt.figure(figsize=(6,5))
    plt.plot(time,voltage,',',color='#1f77b4')
    plt.plot(tVrest,Vrest,color='#1f77b4')
    plt.plot([time[0],time[-1]],[Vrest_avg,Vrest_avg],'--',color='#ff7f0e')
    plt.xlabel('Time [ms]')
    plt.ylabel('Voltage [mV]')
    plt.title('Testing implementation')
    plt.tight_layout()
    plt.show() 
    '''
    
    return Vrest_avg, Vrest_rms


if __name__ == '__main__':
    testmodels = [478513437, 488462965, 478513407] 
    cm         = 1.0
    spikedurat = -40
    iamp       = 0
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
    varyE = [0,13,13.5,14,25,26,26.5,27,28,29,115,116,117,118,119,119.5,120]#[-30,-25,-20,-10.9,-10,-5,0,5,10,11,19.7,20,25,30,40,50]
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
    
    # hmmm...
    changestring =''
    changestring = changestring+'_E'+str(varyE)+'_gdflt'

    ###########################################################
    
    NE = len(varyE)
    zeroind = 3 # Default
    for i in range(NE):
        if varyE[i]==0:
            zeroind = i
    
    Nspikes = numpy.zeros(NE)
    avg_ISI = numpy.zeros(NE)
    rms_ISI = numpy.zeros(NE)
    avg_AP_ampl = numpy.zeros(NE)
    rms_AP_ampl = numpy.zeros(NE)
    avg_AP_mins = numpy.zeros(NE)
    rms_AP_mins = numpy.zeros(NE)
    avg_AP_halfwidth = numpy.zeros(NE)
    rms_AP_halfwidth = numpy.zeros(NE)
    
    outfolder_all = 'figures/Comparemodels/Vrest/'
    plotname_all  = outfolder_all+'%s_Allen_cmfall'%namestringfirst+str(cm)+'_idur%i_'% idur+'_manual_Vrest_vs_Epas.png'
    outfilename_all  = outfolder_all+'%s_Allen_cmfall'%namestringfirst+str(cm)+'_idur%i_'% idur+'_manual_Vrestdiff_vs_Epas.txt'
    print('varyE:',varyE)
    # Set names
    Vrest_avg_all = []
    Vrest_rms_all = []
    for testmodel in testmodels:
        
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
        outfilename = outfolder+'%s_%i_cmfall'%(namestringfirst,testmodel)+str(cm)+'_idur%i_'% idur+'_manual_Vrest_vs_Epas.txt'
        plotname  = outfolder+'%s_%i_cmfall'%(namestringfirst,testmodel)+str(cm)+'_idur%i_'% idur+'_manual_Vrest_vs_Epas.png'
        # make files
        outfile = open(outfilename,'w')
        thr_reached=False
        Vrest_avg_these = numpy.zeros(NE)
        Vrest_rms_these = numpy.zeros(NE)
        for j in range(NE):
            E = varyE[j]
            namestring = namestringfirst+str(E)
            infolder = 'figures/%i/current_idur%i_iamp' % (testmodel,idur) + str(iamp)+'/'
            print('Step ', j+1, ' of', NE)
            filename = infolder+namestring+"_changecmf" + str(cm) + "_everywhere_vinit"+str(v_init)+"_addedRa.txt"
            Vrest_avg_these[j], Vrest_rms_these[j] = manual(filename,idelay,idur,spikedurat)  
        outfile.close()
        Vrest_avg_all.append(Vrest_avg_these)
        Vrest_rms_all.append(Vrest_rms_these)
        
        # Plot results
        plt.figure(figsize=(6,5))
        plt.fill_between(varyE,Vrest_avg_these+Vrest_rms_these,Vrest_avg_these+Vrest_rms_these, alpha=0.3)
        plt.plot(varyE,Vrest_avg_these)
        plt.xlabel(r'$E_{\mathregular{pas}}$ (mV)')
        plt.ylabel(r'$V_{\mathregular{rest}}$ (mV)')
        plt.title(r'$E_{\mathregular{pas}}$ vs $V_{\mathregular{rest}}$, model %i $C_{m,all}$' % testmodel)
        plt.tight_layout()
        plt.savefig(plotname)

Vrest_avg_1 = Vrest_avg_all[0]
Vrest_rms_1 = Vrest_rms_all[0]
Vrest_avg_2 = Vrest_avg_all[1]
Vrest_rms_2 = Vrest_rms_all[1]
Vrest_avg_3 = Vrest_avg_all[2]
Vrest_rms_3 = Vrest_rms_all[2]
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(25,10))

# Model 1
ax1.fill_between(varyE,Vrest_avg_1+Vrest_rms_1,Vrest_avg_1+Vrest_rms_1, alpha=0.3)
ax1.plot(varyE,Vrest_avg_1)
#ax1.axhline(y=0,color='k',linestyle='--',linewidth=1.0)
#ax1.axvline(x=Vrest_avg_1[zeroind],color='k',linestyle='--',linewidth=1.0)
ax1.set_xlabel(r'$E_{\mathregular{pas}}$ (mV)')
ax1.set_ylabel(r'$V_{\mathregular{rest}}$ (mV)')
ax1.set_title(r'$E_{\mathregular{pas}}$ vs $V_{\mathregular{rest}}$, model %i $C_{m,all}$' % testmodels[0])

# Model 2
ax2.fill_between(varyE,Vrest_avg_2+Vrest_rms_2,Vrest_avg_2+Vrest_rms_2, alpha=0.3)
ax2.plot(varyE,Vrest_avg_2)
#ax2.axhline(y=0,color='k',linestyle='--',linewidth=1.0)
#ax2.axvline(x=Vrest_avg_1[zeroind],color='k',linestyle='--',linewidth=1.0)
ax2.set_xlabel(r'$E_{\mathregular{pas}}$ (mV)')
ax2.set_ylabel(r'$V_{\mathregular{rest}}$ (mV)')
ax2.set_title(r'$E_{\mathregular{pas}}$ vs $V_{\mathregular{rest}}$, model %i $C_{m,all}$' % testmodels[1])

# Model 3
ax3.fill_between(varyE,Vrest_avg_3+Vrest_rms_3,Vrest_avg_3+Vrest_rms_3, alpha=0.3)
ax3.plot(varyE,Vrest_avg_3)
#ax3.axhline(y=0,color='k',linestyle='--',linewidth=1.0)
#ax3.axvline(x=Vrest_avg_1[zeroind],color='k',linestyle='--',linewidth=1.0)
ax3.set_xlabel(r'$E_{\mathregular{pas}}$ (mV)')
ax3.set_ylabel(r'$V_{\mathregular{rest}}$ (mV)')
ax3.set_title(r'$E_{\mathregular{pas}}$ vs $V_{\mathregular{rest}}$, model %i $C_{m,all}$' % testmodels[2])

fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(plotname_all)

outfile_all = open(outfilename_all,'w')
outfile_all.write('Epas   Model437   Model965   Model407\n')
Vbase_1 = Vrest_avg_1[zeroind]
Vbase_2 = Vrest_avg_2[zeroind]
Vbase_3 = Vrest_avg_3[zeroind]
for i in range(NE):
    E = varyE[i]
    Vdiff1 = Vrest_avg_1[i]-Vbase_1
    Vdiff2 = Vrest_avg_2[i]-Vbase_2
    Vdiff3 = Vrest_avg_3[i]-Vbase_3
    outfile_all.write('%.2f %.3f %.3f %.3f\n' % (E,Vdiff1,Vdiff2,Vdiff3))
outfile_all.close()
plt.show()