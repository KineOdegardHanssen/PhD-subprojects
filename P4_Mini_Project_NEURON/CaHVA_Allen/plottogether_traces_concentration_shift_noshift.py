import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np

def find_peakind(values):
    N = len(values)
    peakinds = []
    peakvals = []
    for i in range(1,N-1):
        valuesi = values[i]
        if valuesi>values[i-1] and valuesi>values[i+1]:
            peakinds.append(i)
            peakvals.append(valuesi)
    return peakinds, peakvals

def givepeak(t,V):
    N = len(t)
    peaktimes = []
    for i in range(1,N-1):
        if V[i+1]<V[i] and V[i-1]<V[i]:
            peaktimes.append(t[i])
    return peaktimes

def unpack_cahva(filename):
    data = np.loadtxt(filename)
    
    t = data[:, 0]
    v = data[:, 1]
    eca = data[:, 2]
    cai = data[:, 3]
    cao = data[:, 4]
    I_Ca_HVA = data[:, 5]
    g_Ca_HVA = data[:, 6]
    return t, v, eca, cai, cao, I_Ca_HVA, g_Ca_HVA

def unpack_sk(filename):
    data = np.loadtxt(filename)
    
    t = data[:, 0]
    v = data[:, 1]
    eca = data[:, 2]
    cai = data[:, 3]
    cao = data[:, 4]
    I_SK = data[:, 5]
    I_Ca_HVA = data[:, 6]
    g_SK     = data[:, 7]
    g_Ca_HVA = data[:, 8]
    return t, v, eca, cai, cao, I_SK, I_Ca_HVA, g_SK, g_Ca_HVA

############################# No shifts ############################################

def avg_and_rms(x):
    N = len(x)
    avgx = np.mean(x)
    rmsx = 0
    for i in range(N):
        rmsx += (avgx-x[i])**2
    rmsx = np.sqrt(rmsx/(N-1))
    return avgx,rmsx

def Nspikes(time,voltage,idelay,idur,skiptime):
    """Manual approach"""
    
    anfraction = (idur-skiptime)/float(idur)
    tstartan   = idelay+skiptime
    vmax = max(voltage) 
    vmin = min(voltage) 
    deltav = vmax-vmin
    vthr  = -20   # If there is a peak above this value, we count it # Old: vmax-0.15*deltav
    vprev = vthr-40 # A peak never kicks in at initiation, even if I change vthr
    durthr = spikedurat # Height at which we measure the duration. Need to be calibrated by peaks? # Was -20
    Npeaks = 0
    peakmins  = []
    peakvals  = []
    peaktimes = []
    passtimes_up = []
    passvals_up  = []
    passtimes_down = []
    passvals_down  = []
    minalready = False
    for i in range (1,len(voltage)-1):  
        if time[i]<idelay+idur:
            if voltage[i-1]<voltage[i] and voltage[i+1]<voltage[i] and voltage[i]>vthr and time[i]>tstartan:
                peaktimes.append(time[i])
                peakvals.append(voltage[i])
                Npeaks+=1
                minalready = False
            if voltage[i-1]>voltage[i] and voltage[i+1]>voltage[i] and voltage[i]<vthr and time[i]>tstartan and minalready==False and Npeaks>0:
                peakmins.append(voltage[i])
                minalready = True
            if voltage[i]>=durthr and voltage[i-1]<durthr and time[i]>tstartan: # Passing upwards
                tbef = time[i-1]
                taft = time[i]
                Vbef = voltage[i-1]
                Vaft = voltage[i]
                a = (Vaft-Vbef)/(taft-tbef)
                b = Vbef-a*tbef
                tint = (durthr-b)/a
                Vint = a*tint+b
                passtimes_up.append(tint)
                passvals_up.append(Vint) # For plotting
            elif voltage[i]>=durthr and voltage[i+1]<durthr and len(passtimes_up)>0: # Passing downwards
                tbef = time[i]
                taft = time[i+1]
                Vbef = voltage[i]
                Vaft = voltage[i+1]
                a = (Vaft-Vbef)/(taft-tbef)
                b = Vbef-a*tbef
                tint = (durthr-b)/a
                Vint = a*tint+b
                passtimes_down.append(tint)
                passvals_down.append(Vint) # For plotting
        else:
            break
    
    Npeaks /= anfraction # Want frequency
    # Checking if we've got consistent firing:
    if Npeaks!=0:
        if peaktimes[-1]<=(idur/2.+idelay): #Checking if there's no firing in the last half of the stim. interval
            Npeaks=0                        # Is that a proper limit? Last third instead?
    
    # I should probably plot stuff, to make sure... If-test?
    dur = []
    isi = []
    amps = []
    Namps = min([len(peakmins),len(peakvals)])
    Ndur = min([len(passtimes_up),len(passtimes_down)]) # Should be the same
    for i in range(Ndur-1):
        dur.append(passtimes_down[i]-passtimes_up[i])
    for i in range(Namps):
        amps.append(peakvals[i]-peakmins[i])
    for i in range(1,len(peaktimes)):
        isi.append(peaktimes[i]-peaktimes[i-1])
    #isi      = isi[:-1]
    #peakvals = peakvals[:-1]
    #peakmins = peakmins[:-1]
    time_peakvals = peaktimes#[:-1]
    
    '''
    plt.figure(figsize=(6,5))
    plt.plot(time,voltage,',')
    plt.plot(time_peakvals,peakvals,'o',label='peaks')
    plt.plot(passtimes_up,passvals_up,'o',label='dur basis, up')
    plt.plot(passtimes_down,passvals_down,'o',label='dur basis, down')
    plt.xlabel('Time [ms]')
    plt.ylabel('Voltage [mV]')
    plt.title('Testing implementation')
    plt.legend(loc='lower left')
    plt.tight_layout()
    plt.show() 
    '''
    
    ## Avg and rms:
    amps_avg, amps_rms = avg_and_rms(amps)
    peakmins_avg, peakmins_rms = avg_and_rms(peakmins)
    peakvals_avg, peakvals_rms = avg_and_rms(peakvals)
    dur_avg, dur_rms = avg_and_rms(dur)
    isi_avg, isi_rms = avg_and_rms(isi)
    
    return Npeaks #, peaktimes, peakmins_avg, peakmins_rms, peakvals_avg,  peakvals_rms, dur_avg, dur_rms, isi_avg, isi_rms, isi, dur, amps_avg, amps_rms

##
def heightdiff(time,concentration,idelay,idur,skiptime):
    """Manual approach"""
    
    anfraction = (idur-skiptime)/float(idur)
    tstartan   = idelay+skiptime
    cmax = max(concentration) 
    cmin = min(concentration) 
    vthr  = -20   # If there is a peak above this value, we count it
    durthr = spikedurat # Height at which we measure the duration. Need to be calibrated by peaks? # Was -20
    Npeaks = 0
    peakmins  = []
    peakvals  = []
    peaktimes = []
    passtimes_up = []
    passvals_up  = []
    passtimes_down = []
    passvals_down  = []
    minalready = False
    for i in range (1,len(concentration)-1):  
        if time[i]<idelay+idur:
            if concentration[i-1]<concentration[i] and concentration[i+1]<concentration[i] and time[i]>tstartan:
                peaktimes.append(time[i])
                peakvals.append(concentration[i])
                Npeaks+=1
                minalready = False
            if concentration[i-1]>concentration[i] and concentration[i+1]>concentration[i] and time[i]>tstartan and minalready==False and Npeaks>0:
                peakmins.append(concentration[i])
                minalready = True
            if concentration[i]>=durthr and concentration[i-1]<durthr and time[i]>tstartan: # Passing upwards
                tbef = time[i-1]
                taft = time[i]
                Vbef = concentration[i-1]
                Vaft = concentration[i]
                a = (Vaft-Vbef)/(taft-tbef)
                b = Vbef-a*tbef
                tint = (durthr-b)/a
                Vint = a*tint+b
                passtimes_up.append(tint)
                passvals_up.append(Vint) # For plotting
            elif concentration[i]>=durthr and concentration[i+1]<durthr and len(passtimes_up)>0: # Passing downwards
                tbef = time[i]
                taft = time[i+1]
                Vbef = concentration[i]
                Vaft = concentration[i+1]
                a = (Vaft-Vbef)/(taft-tbef)
                b = Vbef-a*tbef
                tint = (durthr-b)/a
                Vint = a*tint+b
                passtimes_down.append(tint)
                passvals_down.append(Vint) # For plotting
        else:
            break
    
    Npeaks /= anfraction # Want frequency
    # Checking if we've got consistent firing:
    if Npeaks!=0:
        if peaktimes[-1]<=(idur/2.+idelay): #Checking if there's no firing in the last half of the stim. interval
            Npeaks=0                        # Is that a proper limit? Last third instead?
    
    # I should probably plot stuff, to make sure... If-test?
    dur = []
    isi = []
    amps = []
    Namps = min([len(peakmins),len(peakvals)])
    Ndur = min([len(passtimes_up),len(passtimes_down)]) # Should be the same
    for i in range(Ndur-1):
        dur.append(passtimes_down[i]-passtimes_up[i])
    for i in range(Namps):
        amps.append(peakvals[i]-peakmins[i])
    for i in range(1,len(peaktimes)):
        isi.append(peaktimes[i]-peaktimes[i-1])
    #isi      = isi[:-1]
    #peakvals = peakvals[:-1]
    #peakmins = peakmins[:-1]
    time_peakvals = peaktimes#[:-1]
    
    '''
    plt.figure(figsize=(6,5))
    plt.plot(time,voltage,',')
    plt.plot(time_peakvals,peakvals,'o',label='peaks')
    plt.plot(passtimes_up,passvals_up,'o',label='dur basis, up')
    plt.plot(passtimes_down,passvals_down,'o',label='dur basis, down')
    plt.xlabel('Time [ms]')
    plt.ylabel('Voltage [mV]')
    plt.title('Testing implementation')
    plt.legend(loc='lower left')
    plt.tight_layout()
    plt.show() 
    '''
    
    ## Avg and rms:
    amps_avg, amps_rms = avg_and_rms(amps)
    peakmins_avg, peakmins_rms = avg_and_rms(peakmins)
    peakvals_avg, peakvals_rms = avg_and_rms(peakvals)
    dur_avg, dur_rms = avg_and_rms(dur)
    isi_avg, isi_rms = avg_and_rms(isi)
    
    percentdiff = 100*(peakvals_avg-peakmins_avg)/peakmins_avg
    
    return amps_avg, amps_rms, percentdiff, peakmins_avg, peakmins_rms, peakvals_avg,  peakvals_rms #, peaktimes, peakmins_avg, peakmins_rms, peakvals_avg,  peakvals_rms, dur_avg, dur_rms, isi_avg, isi_rms, isi, dur, amps_avg, amps_rms


idur = 1000
iamp = 0.02 # 0.04 # 

plot_all  = False # True # 

gbarbase = 0.00014931667074610222
gcahvas = [0.2]
gcahvas_label = [1.0]
Ng = len(gcahvas)

colors_gCa = []
for i in range(Ng):
    colors_gCa.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))

folderbase    = 'Results/Soma10/current_idur'+str(idur)+'_iamp'+str(iamp)+'/'

ts = []
Vs = []
ecas = []
cais = []
caos = []
I_Ca_HVAs = []
g_Ca_HVAs = []

for i in range(Ng):
    gcahva_this = gcahvas[i]
    filename = folderbase+'_gCaHVA'+str(gcahva_this)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V_cainfo.txt'
    
    print('filename:',filename)
    t, v, eca, cai, cao, I_Ca_HVA, g_Ca_HVA = unpack_cahva(filename)
    g_Ca_HVA = g_Ca_HVA/float(gcahva_this*gbarbase)
    N = len(t)
    ts.append(t)
    Vs.append(v)
    ecas.append(eca)
    cais.append(cai)
    caos.append(cao)
    I_Ca_HVAs.append(I_Ca_HVA)
    g_Ca_HVAs.append(g_Ca_HVA)

t1 = ts[0]
V1 = Vs[0]
eca1 = ecas[0]
cai1 = cais[0]
cao1 = caos[0]
I_Ca_HVA1 = I_Ca_HVAs[0]
g_Ca_HVA1 = g_Ca_HVAs[0]


######################## SK ####################################

#iamp     = 0.04
gsk      = 1.0
gcahva   = 0.2
gsk_base = 0.0028175455472127641

gsks = [1,2,3,5,10]

folderbase    = 'Results/Soma10/current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
filename_sk1  = 'SK_Allen/'+folderbase+'_gSK'+str(gsk)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V_eca.txt'
print('filename:',filename_sk1)


t_sk1, V_sk1, eca_sk1, cai_sk1, cao_sk1, I_SK_sk1, I_Ca_HVA_sk1, g_SK_sk1, g_Ca_HVA_sk1 = unpack_sk(filename_sk1)

g_SK_sk1 = g_SK_sk1/float(gsk*gsk_base)


gcahva_base = 0.00014931667074610222
g_Ca_HVA_sk1  = g_Ca_HVA_sk1/float(gcahva*gcahva_base)

Ng = 6
colors_gSK = []
for i in range(Ng):
    colors_gSK.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))

#'''
idelay   = 100
skiptime = 500
spikedurat = -40
Nspikes1 = Nspikes(t_sk1, V_sk1,idelay,idur,skiptime)
print('Nspikes, 1gSK:',Nspikes1)

amps_avg_gsk1, amps_rms_gsk1, percentdiff1, peakmins_avg1, peakmins_rms1, peakvals_avg1,  peakvals_rms1 = heightdiff(t_sk1,cai_sk1,idelay,idur,skiptime)


####################################################################################


idur = 1000
iamp = 0.02

plotitall = False # True # 

gcahva = 0.2
gcahva_label = 1.0
gcahva_base  = 0.00014931667074610222
Vshifts = [0]
NV = len(Vshifts)


colors = []
for i in range(NV):
    colors.append((1.0-i/float(NV),0,i/float(NV),1.0-abs(0.5-i/float(NV))))

ts = []
Vs = []
ecas = []
cais = []
caos = []
I_Ca_HVAs = []
g_Ca_HVAs = []

for i in range(NV):
    folder = 'Shift_gCaHVA_V/Vshift_'+str(Vshifts[i])+'/Results/Soma10/current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    filename = folder+'_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V_cainfo.txt'
    print('filename:',filename)
    
    t, v, eca, cai, cao, I_Ca_HVA, g_Ca_HVA = unpack_cahva(filename)
    g_Ca_HVA = g_Ca_HVA/float(gcahva*gcahva_base)
    N = len(t)
    #ax1.plot(t[int(N/2):], V[int(N/2):],label=r'%s$\bar{g}$' % gcahvas_label[i])
    ts.append(t)
    Vs.append(v)
    ecas.append(eca)
    cais.append(cai)
    caos.append(cao)
    I_Ca_HVAs.append(I_Ca_HVA)
    g_Ca_HVAs.append(g_Ca_HVA)
#ax1.legend(loc='upper left',ncol=1,fontsize=11)

t1_shift = ts[0]
V1_shift = Vs[0]
eca1_shift = ecas[0]
cai1_shift = cais[0]
cao1_shift = caos[0]
I_Ca_HVA1_shift = I_Ca_HVAs[0]
g_Ca_HVA1_shift = g_Ca_HVAs[0]


######################## SK ####################################
gsk      = 1.0

ts = []
Vs = []
ecas = []
cais = []
caos = []
I_SKs = []
I_Ca_HVAs = []
g_SKs     = []
g_Ca_HVAs = []
for i in range(NV):
    folder = 'SK_Allen/Shift_gCaHVA_V/Vshift_'+str(Vshifts[i])+'/Results/Soma10/current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    filename = folder+'_gSK'+str(gsk)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V_writecainfo.txt'
    print('filename:',filename)
    t, V, eca_sk, cai_sk, cao_sk, I_SK_sk, I_Ca_HVA_sk, g_SK_sk, g_Ca_HVA_sk = unpack_sk(filename)
    N = len(t)
    #ax1.plot(t[int(N/2):], V[int(N/2):],label=r'%s$\bar{g}$' % gcahvas_label[i])
    ts.append(t)
    Vs.append(V)
    ecas.append(eca_sk)
    cais.append(cai_sk)
    caos.append(cao_sk)
    I_SKs.append(I_SK_sk)
    g_SKs.append(g_SK_sk)
    I_Ca_HVAs.append(I_Ca_HVA_sk)
    g_Ca_HVAs.append(g_Ca_HVA_sk)

t_sk1_shift = ts[0]
V_sk1_shift = Vs[0]
eca_sk1_shift = ecas[0]
cai_sk1_shift = cais[0]
cao_sk1_shift = caos[0]
I_SK_sk1_shift = I_SKs[0]
g_SK_sk1_shift = g_SKs[0]
I_Ca_HVA_sk1_shift = I_Ca_HVAs[0]
g_Ca_HVA_sk1_shift = g_Ca_HVAs[0]

gsk_base = 0.0028175455472127641
g_SK_sk1_shift = g_SK_sk1_shift/float(gsk*gsk_base)
g_Ca_HVA_sk1_shift = g_Ca_HVA_sk1_shift/float(gcahva*gcahva_base)


################## Plotting total tracde and [Ca]in: ##################
fig = plt.figure(figsize=(10,8))
    
gs = gridspec.GridSpec(2, 4)
    
ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[1, 0:2])
ax4 = plt.subplot(gs[1, 2:4])

ax1.set_title('A',loc='left')
ax2.set_title('B',loc='left')
ax3.set_title('C',loc='left')
ax4.set_title('D',loc='left')   


t1diff = t1-t1_shift
V1diff = V1-V1_shift
tsk1diff = t_sk1-t_sk1_shift
Vsk1diff = V_sk1-V_sk1_shift

cai1diff = cai1-cai1_shift
caisk1diff = cai_sk1-cai_sk1_shift

ax1.plot(t1, t1diff)
ax1.set_xlabel(r'$t$ (ms)')
ax1.set_ylabel(r'$V$ (mV)')
ax1.set_title('Trace, CaHVA only')
#ax1.set_ylim(bottom=-85)

ax2.plot(t1, cai1diff)
ax2.set_xlabel(r'$t$ (ms)')
ax2.set_ylabel(r'$[\mathregular{Ca}^{2+}]_\mathregular{in}$ (mM)')
ax2.set_title(r'$[\mathregular{Ca}^{2+}]_\mathregular{in}$, CaHVA only')

ax3.plot(t_sk1, Vsk1diff)
ax3.set_xlabel(r'$t$ (ms)')
ax3.set_ylabel(r'$V$ (mV)')
ax3.set_title('Trace, CaHVA and SK')
#ax3.set_ylim(bottom=-90)

ax4.plot(t_sk1, caisk1diff)
ax4.set_xlabel(r'$t$ (ms)')
ax4.set_ylabel(r'$[\mathregular{Ca}^{2+}]_\mathregular{in}$ (mM)')
ax4.set_title(r'$[\mathregular{Ca}^{2+}]_\mathregular{in}$, CaHVA and SK')

plt.tight_layout()
plt.savefig('Results/Soma10/Compare/trace_concentration_idur'+str(idur)+'_iamp'+str(iamp)+'_Ca_SK_SHIFTS_both_DIFF.png')

