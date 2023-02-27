import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np

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

idur = 1000
iamp = 0.02

plotitall = False # True # 

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

fig = plt.figure(figsize=(10,15),dpi=300)#(figsize=(8,3),dpi=300)

gs = gridspec.GridSpec(3, 4)
    
ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])    
ax3 = plt.subplot(gs[1, 0:2])
ax4 = plt.subplot(gs[1, 2:4])    
ax5 = plt.subplot(gs[2, 0:2])
ax6 = plt.subplot(gs[2, 2:4])
    
fig.suptitle(r'$I=$%s nA' % str(iamp),fontsize=20)

ax1.set_title(r'A',loc='left',fontsize=16)
ax2.set_title(r'B',loc='left',fontsize=16)
ax3.set_title(r'C',loc='left',fontsize=16)
ax4.set_title(r'D',loc='left',fontsize=16)
ax5.set_title(r'E',loc='left',fontsize=16)
ax6.set_title(r'F',loc='left',fontsize=16)
ax1.set_title(r'Voltage trace',fontsize=16)#loc='right',
ax2.set_title(r'$E_\mathregular{Ca}$',fontsize=16)#loc='right',
ax3.set_title(r'Voltage trace difference',fontsize=16)#loc='right',
ax4.set_title(r'$\left[\mathregular{Ca}^{2+}\right]_\mathregular{in}$',fontsize=16)#loc='right',
ax5.set_title(r'$I_\mathregular{Ca}$',fontsize=16)#loc='right',
ax6.set_title(r'$g_\mathregular{Ca}$',fontsize=16)#loc='right',
mycolors = ['tab:orange','tab:green','tab:red','tab:purple','tab:brown']
    
ax1.set_xlabel(r'$t$ (ms)',fontsize=14)
ax1.set_ylabel(r'$V$ (mV)',fontsize=14)
ax2.set_xlabel(r'$t$ (ms)',fontsize=14)
ax2.set_ylabel(r'$V$ (mV)',fontsize=14)
ax3.set_xlabel(r'$t$ (ms)',fontsize=14)
ax3.set_ylabel(r'$\Delta V$ (mV)',fontsize=14)
ax4.set_xlabel(r'$t$ (ms)',fontsize=14)
ax4.set_ylabel(r'$\left[\mathregular{Ca}^{2+}\right]_\mathregular{in}$ (mM)',fontsize=14)
ax5.set_xlabel(r'$t$ (ms)',fontsize=14)
ax5.set_ylabel(r'$I_\mathregular{Ca}$ (nA)',fontsize=14)
ax6.set_xlabel(r'$t$ (ms)',fontsize=14)
ax6.set_ylabel(r'$g_\mathregular{Ca}/\bar{g}_\mathregular{Ca}$',fontsize=14)

gbarbase = 0.00014931667074610222
gcahvas = [0.2,0.5,1.0,1.5]#,2.0]
gcahvas_label = [1.0,2.5,5.0,7.5]#,10.0]
Ng = len(gcahvas)

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
    
    t, v, eca, cai, cao, I_Ca_HVA, g_Ca_HVA = unpack_cahva(filename)
    g_Ca_HVA = g_Ca_HVA/float(gcahva_this*gbarbase)
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

t1 = ts[0] # 0.2
t2 = ts[1] # 0.5
t3 = ts[2] # 1.0
t4 = ts[3] # 1.5
#t5 = ts[4] # 2.0
V1 = Vs[0]
V2 = Vs[1]
V3 = Vs[2]
V4 = Vs[3]
#V5 = Vs[4]
eca1 = ecas[0] # 0.2
eca2 = ecas[1] # 0.5
eca3 = ecas[2] # 1.0
eca4 = ecas[3] # 1.5
###
cai1 = cais[0] # 0.2
cai2 = cais[1] # 0.5
cai3 = cais[2] # 1.0
cai4 = cais[3] # 1.5
###
cao1 = caos[0] # 0.2
cao2 = caos[1] # 0.5
cao3 = caos[2] # 1.0
cao4 = caos[3] # 1.5
###
I_Ca_HVA1 = I_Ca_HVAs[0] # 0.2
I_Ca_HVA2 = I_Ca_HVAs[1] # 0.5
I_Ca_HVA3 = I_Ca_HVAs[2] # 1.0
I_Ca_HVA4 = I_Ca_HVAs[3] # 1.5
###
g_Ca_HVA1 = g_Ca_HVAs[0] # 0.2
g_Ca_HVA2 = g_Ca_HVAs[1] # 0.5
g_Ca_HVA3 = g_Ca_HVAs[2] # 1.0
g_Ca_HVA4 = g_Ca_HVAs[3] # 1.5


coshift = -300
ishift1 = 1900
ishift2 = 523
ishift3 = -516+coshift
ishift4 = 5573+coshift

istart = int(5*N/8)+5549
iend   = int(3*N/4)-11458

istart1 = istart-ishift1
iend1   = iend-ishift1
istart2 = istart-ishift2
iend2   = iend-ishift2
istart3 = istart-ishift3
iend3   = iend-ishift3
istart4 = istart-ishift4
iend4   = iend-ishift4
#istart5 = istart+3874#
#iend5   = iend-3283
#ishift5 = 615

t1_shift = t1[istart1:iend1]
V1_shifted_CaHVA = V1[istart1:iend1]
V2_shifted_CaHVA = V2[istart2:iend2]
V3_shifted_CaHVA = V3[istart3:iend3]
V4_shifted_CaHVA = V4[istart4:iend4]
#V5_shifted_CaHVA = V5[istart5:iend5]
eca1_shifted_CaHVA = eca1[istart1:iend1] # 0.2
eca2_shifted_CaHVA = eca2[istart2:iend2] # 0.5
eca3_shifted_CaHVA = eca3[istart3:iend3] # 1.0
eca4_shifted_CaHVA = eca4[istart4:iend4] # 1.5
###
cai1_shifted_CaHVA = cai1[istart1:iend1] # 0.2
cai2_shifted_CaHVA = cai2[istart2:iend2] # 0.5
cai3_shifted_CaHVA = cai3[istart3:iend3] # 1.0
cai4_shifted_CaHVA = cai4[istart4:iend4] # 1.5
###
I_Ca_HVA1_shifted_CaHVA = I_Ca_HVA1[istart1:iend1] # 0.2
I_Ca_HVA2_shifted_CaHVA = I_Ca_HVA2[istart2:iend2] # 0.5
I_Ca_HVA3_shifted_CaHVA = I_Ca_HVA3[istart3:iend3] # 1.0
I_Ca_HVA4_shifted_CaHVA = I_Ca_HVA4[istart4:iend4] # 1.5
###
g_Ca_HVA1_shifted_CaHVA = g_Ca_HVA1[istart1:iend1] # 0.2
g_Ca_HVA2_shifted_CaHVA = g_Ca_HVA2[istart2:iend2] # 0.5
g_Ca_HVA3_shifted_CaHVA = g_Ca_HVA3[istart3:iend3] # 1.0
g_Ca_HVA4_shifted_CaHVA = g_Ca_HVA4[istart4:iend4] # 1.5

deltaV1 = V2_shifted_CaHVA-V1_shifted_CaHVA
deltaV2 = V3_shifted_CaHVA-V1_shifted_CaHVA
deltaV3 = V4_shifted_CaHVA-V1_shifted_CaHVA
#deltaV4 = V5_shifted_CaHVA-V1_shifted_CaHVA

print('len(V1_shifted_CaHVA):',len(V1_shifted_CaHVA))
print('len(V2_shifted_CaHVA):',len(V2_shifted_CaHVA))
print('len(V3_shifted_CaHVA):',len(V3_shifted_CaHVA))
print('len(V4_shifted_CaHVA):',len(V4_shifted_CaHVA))
#print('len(V5_shifted_CaHVA):',len(V5_shifted_CaHVA))


#'''
ax1.plot(t1[istart1+ishift1:iend1+ishift1], V1_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[0]))
ax1.plot(t2[istart2+ishift2:iend2+ishift2], V2_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[1]))
ax1.plot(t3[istart3+ishift3:iend3+ishift3], V3_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[2]))
ax1.plot(t4[istart4+ishift4:iend4+ishift4], V4_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[3]))
#ax1.plot(t5[istart5-ishift5:iend5-ishift5], V5[istart5:iend5],label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[4]))
ax1.set_xlim(left=t1[istart1+ishift1],right=t1[iend1+ishift1])
#ax1.axis([723,746,-78.2,51.1])
ax1.legend(loc='upper left',ncol=1,fontsize=11)

ax2.plot(t1[istart1+ishift1:iend1+ishift1], eca1_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[0]))
ax2.plot(t2[istart2+ishift2:iend2+ishift2], eca2_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[1]))
ax2.plot(t3[istart3+ishift3:iend3+ishift3], eca3_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[2]))
ax2.plot(t4[istart4+ishift4:iend4+ishift4], eca4_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[3]))
#ax1.plot(t5[istart5-ishift5:iend5-ishift5], V5[istart5:iend5],label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[4]))
ax2.set_xlim(left=t1[istart1+ishift1],right=t1[iend1+ishift1])
#ax1.axis([723,746,-78.2,51.1])
ax2.legend(loc='upper left',ncol=1,fontsize=11)

ax3.plot(t1[istart1+ishift1:iend1+ishift1], deltaV1,label=r'$V_{%s\bar{g}_\mathregular{CaHVA}}$-$V_{%s\bar{g}_\mathregular{CaHVA}}$' % (str(gcahvas_label[1]),str(gcahvas_label[0])),color=mycolors[0])
ax3.plot(t1[istart1+ishift1:iend1+ishift1], deltaV2,label=r'$V_{%s\bar{g}_\mathregular{CaHVA}}$-$V_{%s\bar{g}_\mathregular{CaHVA}}$' % (str(gcahvas_label[2]),str(gcahvas_label[0])),color=mycolors[1])
ax3.plot(t1[istart1+ishift1:iend1+ishift1], deltaV3,label=r'$V_{%s\bar{g}_\mathregular{CaHVA}}$-$V_{%s\bar{g}_\mathregular{CaHVA}}$' % (str(gcahvas_label[3]),str(gcahvas_label[0])),color=mycolors[2])
#ax3.plot(t1_shift, deltaV4,label=r'$V_{%s\bar{g}_\mathregular{CaHVA}}$-$V_{%s\bar{g}_\mathregular{CaHVA}}$' % (str(gcahvas_label[4]),str(gcahvas_label[0])),color=mycolors[3])
ax3.set_xlim(left=t1[istart1+ishift1],right=t1[iend1+ishift1])
#ax3.axis([720.109,749.359,-78.2,51.1])
ax3.legend(loc='lower right',ncol=1,fontsize=11)

ax4.plot(t1[istart1+ishift1:iend1+ishift1], cai1_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[0]))
ax4.plot(t2[istart2+ishift2:iend2+ishift2], cai2_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[1]))
ax4.plot(t3[istart3+ishift3:iend3+ishift3], cai3_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[2]))
ax4.plot(t4[istart4+ishift4:iend4+ishift4], cai4_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[3]))
#ax1.plot(t5[istart5-ishift5:iend5-ishift5], V5[istart5:iend5],label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[4]))
ax4.set_xlim(left=t1[istart1+ishift1],right=t1[iend1+ishift1])
#ax1.axis([723,746,-78.2,51.1])
ax4.legend(loc='upper left',ncol=1,fontsize=11)

ax5.plot(t1[istart1+ishift1:iend1+ishift1], I_Ca_HVA1_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[0]))
ax5.plot(t2[istart2+ishift2:iend2+ishift2], I_Ca_HVA2_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[1]))
ax5.plot(t3[istart3+ishift3:iend3+ishift3], I_Ca_HVA3_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[2]))
ax5.plot(t4[istart4+ishift4:iend4+ishift4], I_Ca_HVA4_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[3]))
#ax1.plot(t5[istart5-ishift5:iend5-ishift5], V5[istart5:iend5],label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[4]))
ax5.set_xlim(left=t1[istart1+ishift1],right=t1[iend1+ishift1])
#ax1.axis([723,746,-78.2,51.1])
ax5.legend(loc='lower left',ncol=1,fontsize=11)

ax6.plot(t1[istart1+ishift1:iend1+ishift1], g_Ca_HVA1_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[0]))
ax6.plot(t2[istart2+ishift2:iend2+ishift2], g_Ca_HVA2_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[1]))
ax6.plot(t3[istart3+ishift3:iend3+ishift3], g_Ca_HVA3_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[2]))
ax6.plot(t4[istart4+ishift4:iend4+ishift4], g_Ca_HVA4_shifted_CaHVA,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[3]))
#ax1.plot(t5[istart5-ishift5:iend5-ishift5], V5[istart5:iend5],label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[4]))
ax6.set_xlim(left=t1[istart1+ishift1],right=t1[iend1+ishift1])
#ax1.axis([723,746,-78.2,51.1])
ax6.legend(loc='upper left',ncol=1,fontsize=11)

fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig('Results/Soma10/Compare/trace_cainfo_idur'+str(idur)+'_iamp'+str(iamp)+'_CaHVA_aligned.png')

''' # Eca and cai does indeed vary for base+CaHVA
plt.figure()
plt.plot(t1,eca1)
plt.plot(t2,eca2)
plt.plot(t3,eca3)
plt.plot(t4,eca4)
plt.title('ECa')
plt.show()


plt.figure()
plt.plot(t1,cai1)
plt.plot(t2,cai2)
plt.plot(t3,cai3)
plt.plot(t4,cai4)
plt.title('Cai')
plt.show()
'''
###

######################## SK ####################################

iamp     = 0.04
gsk      = 1.0
gsk2     = 2.0
gsk3     = 3.0
gsk4     = 5.0
gsk5     = 10.0
gcahva   = 0.2
gcahva2  = 1.0
gcahva_plot  = 1.0
gcahva_plot2 = 5.0
gsk_base = 0.0028175455472127641

folderbase    = 'Results/Soma10/current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
filename_nosk = 'SK_Allen/'+folderbase+'_gSK0.0p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V_eca.txt'
filename_sk1  = 'SK_Allen/'+folderbase+'_gSK'+str(gsk)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V_eca.txt'
filename_sk2  = 'SK_Allen/'+folderbase+'_gSK'+str(gsk2)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V_eca.txt'
filename_sk3  = 'SK_Allen/'+folderbase+'_gSK'+str(gsk3)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V_eca.txt'
filename_sk4  = 'SK_Allen/'+folderbase+'_gSK'+str(gsk4)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V_eca.txt'
filename_sk5  = 'SK_Allen/'+folderbase+'_gSK'+str(gsk5)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V_eca.txt'


t_sk1, V_sk1, eca_sk1, cai_sk1, cao_sk1, I_SK_sk1, I_Ca_HVA_sk1, g_SK_sk1, g_Ca_HVA_sk1 = unpack_sk(filename_sk1)
t_sk2, V_sk2, eca_sk2, cai_sk2, cao_sk2, I_SK_sk2, I_Ca_HVA_sk2, g_SK_sk2, g_Ca_HVA_sk2 = unpack_sk(filename_sk2)
t_sk3, V_sk3, eca_sk3, cai_sk3, cao_sk3, I_SK_sk3, I_Ca_HVA_sk3, g_SK_sk3, g_Ca_HVA_sk3 = unpack_sk(filename_sk3)
t_sk4, V_sk4, eca_sk4, cai_sk4, cao_sk4, I_SK_sk4, I_Ca_HVA_sk4, g_SK_sk4, g_Ca_HVA_sk4 = unpack_sk(filename_sk4)
t_sk5, V_sk5, eca_sk5, cai_sk5, cao_sk5, I_SK_sk5, I_Ca_HVA_sk5, g_SK_sk5, g_Ca_HVA_sk5 = unpack_sk(filename_sk5)
t_nosk, V_nosk, eca_nosk, cai_nosk, cao_nosk, I_SK_nosk, I_Ca_HVA_nosk, g_SK_nosk, g_Ca_HVA_nosk = unpack_sk(filename_nosk)

g_SK_sk1 = g_SK_sk1/float(gsk*gsk_base)
g_SK_sk2 = g_SK_sk2/float(gsk2*gsk_base)
g_SK_sk3 = g_SK_sk3/float(gsk3*gsk_base)
g_SK_sk4 = g_SK_sk4/float(gsk4*gsk_base)
g_SK_sk5 = g_SK_sk5/float(gsk5*gsk_base)

#'''
idelay   = 100
skiptime = 500
spikedurat = -40
Nspikes0 = Nspikes(t_nosk, V_nosk,idelay,idur,skiptime)
Nspikes1 = Nspikes(t_sk1, V_sk1,idelay,idur,skiptime)
Nspikes2 = Nspikes(t_sk2, V_sk2,idelay,idur,skiptime)
Nspikes5 = Nspikes(t_sk5, V_sk5,idelay,idur,skiptime)

print('Nspikes, 0gSK:',Nspikes0)
print('Nspikes, 1gSK:',Nspikes1)
print('Nspikes, 2gSK:',Nspikes2)
print('Nspikes, 10gSK:',Nspikes5)

amps_avg_gsk0, amps_rms_gsk0, percentdiff0, peakmins_avg0, peakmins_rms0, peakvals_avg0,  peakvals_rms0 = heightdiff(t_nosk,cai_nosk,idelay,idur,skiptime)
amps_avg_gsk1, amps_rms_gsk1, percentdiff1, peakmins_avg1, peakmins_rms1, peakvals_avg1,  peakvals_rms1 = heightdiff(t_sk1,cai_sk1,idelay,idur,skiptime)
amps_avg_gsk2, amps_rms_gsk2, percentdiff2, peakmins_avg2, peakmins_rms2, peakvals_avg2,  peakvals_rms2 = heightdiff(t_sk2,cai_sk2,idelay,idur,skiptime)
amps_avg_gsk5, amps_rms_gsk5, percentdiff5, peakmins_avg5, peakmins_rms5, peakvals_avg5,  peakvals_rms5 = heightdiff(t_sk5,cai_sk5,idelay,idur,skiptime)

'''
print('amps, 0gSK:', amps_avg_gsk0, '+/-', amps_rms_gsk0)
print('amps, 1gSK:', amps_avg_gsk1, '+/-', amps_rms_gsk1)
print('amps, 2gSK:', amps_avg_gsk2, '+/-', amps_rms_gsk2)
print('amps, 10gSK:', amps_avg_gsk5, '+/-', amps_rms_gsk5)

print('mins, 0gSK:', peakmins_avg0, '+/-', peakmins_rms0)
print('mins, 1gSK:', peakmins_avg1, '+/-', peakmins_rms1)
print('mins, 2gSK:', peakmins_avg2, '+/-', peakmins_rms2)
print('mins, 10gSK:', peakmins_avg5, '+/-', peakmins_rms5)

print('maxes, 0gSK:', peakvals_avg0, '+/-', peakvals_rms0)
print('maxes, 1gSK:', peakvals_avg1, '+/-', peakvals_rms1)
print('maxes, 2gSK:', peakvals_avg2, '+/-', peakvals_rms2)
print('maxes, 10gSK:', peakvals_avg5, '+/-', peakvals_rms5)
'''

print('% diff, 0gSK:', percentdiff0)
print('% diff, 1gSK:', percentdiff1)
print('% diff, 2gSK:', percentdiff2)
print('% diff, 10gSK:', percentdiff5)


fig = plt.figure(figsize=(10,5))#,dpi=300)#(figsize=(8,3),dpi=300)

gs = gridspec.GridSpec(1, 4)
    
ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])   


ax1.plot(t_nosk, V_nosk, label=r'0.0$\bar{g}_\mathregular{SK}$')
ax1.plot(t_sk1, V_sk1, label=r'%.1f$\bar{g}_\mathregular{SK}$' % gsk)
ax1.plot(t_sk2, V_sk2, label=r'%.1f$\bar{g}_\mathregular{SK}$' % gsk2)
ax1.plot(t_sk5, V_sk5, label=r'%.1f$\bar{g}_\mathregular{SK}$' % gsk5)
ax1.set_xlabel(r'$t$ (ms)')
ax1.set_ylabel(r'$V$ (mV)')
ax1.legend(loc='lower center',ncol=4,fontsize=7)

ax2.plot(t_nosk, cai_nosk, label=r'0.0$\bar{g}_\mathregular{SK}$')
ax2.plot(t_sk1, cai_sk1, label=r'%.1f$\bar{g}_\mathregular{SK}$' % gsk)
ax2.plot(t_sk2, cai_sk2, label=r'%.1f$\bar{g}_\mathregular{SK}$' % gsk2)
ax2.plot(t_sk5, cai_sk5, label=r'%.1f$\bar{g}_\mathregular{SK}$' % gsk5)
ax2.set_xlabel(r'$t$ (ms)')
ax2.set_ylabel(r'$[\mathregular{Ca}^{2+}]_\mathregular{in}$ (mM)')
ax2.legend(loc='upper right',ncol=2)

plt.tight_layout()
plt.savefig('Results/Soma10/Compare/trace_concentration_idur'+str(idur)+'_iamp'+str(iamp)+'_SK.png')
#plt.show()
#'''
fig = plt.figure(figsize=(12,20),dpi=300)#(figsize=(8,3),dpi=300)

gs = gridspec.GridSpec(4, 4)
    
ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])    
ax3 = plt.subplot(gs[1, 0:2])
ax4 = plt.subplot(gs[1, 2:4])    
ax5 = plt.subplot(gs[2, 0:2])
ax6 = plt.subplot(gs[2, 2:4])    
ax7 = plt.subplot(gs[3, 0:2])
ax8 = plt.subplot(gs[3, 2:4])
    
fig.suptitle(r'$I=$%s nA' % str(iamp),fontsize=20)

ax1.set_title(r'A',loc='left',fontsize=16)
ax2.set_title(r'B',loc='left',fontsize=16)
ax3.set_title(r'C',loc='left',fontsize=16)
ax4.set_title(r'D',loc='left',fontsize=16)
ax5.set_title(r'E',loc='left',fontsize=16)
ax6.set_title(r'F',loc='left',fontsize=16)
ax7.set_title(r'G',loc='left',fontsize=16)
ax8.set_title(r'H',loc='left',fontsize=16)
ax1.set_title(r'Voltage trace',fontsize=16)#loc='right',
ax2.set_title(r'$E_\mathregular{Ca}$',fontsize=16)#loc='right',
ax3.set_title(r'Voltage trace difference',fontsize=16)#loc='right',
ax4.set_title(r'$\left[\mathregular{Ca}^{2+}\right]_\mathregular{in}$',fontsize=16)#loc='right',
ax5.set_title(r'$I_\mathregular{Ca}$',fontsize=16)#loc='right',
ax6.set_title(r'$g_\mathregular{Ca}$',fontsize=16)#loc='right',
ax7.set_title(r'$I_\mathregular{SK}$',fontsize=16)#loc='right',
ax8.set_title(r'$g_\mathregular{SK}/\bar{g}_\mathregular{SK}$',fontsize=16)#loc='right',
mycolors = ['tab:orange','tab:green','tab:red','tab:purple','tab:brown']
    
ax1.set_xlabel(r'$t$ (ms)',fontsize=14)
ax1.set_ylabel(r'$V$ (mV)',fontsize=14)
ax2.set_xlabel(r'$t$ (ms)',fontsize=14)
ax2.set_ylabel(r'$V$ (mV)',fontsize=14)
ax3.set_xlabel(r'$t$ (ms)',fontsize=14)
ax3.set_ylabel(r'$\Delta V$ (mV)',fontsize=14)
ax4.set_xlabel(r'$t$ (ms)',fontsize=14)
ax4.set_ylabel(r'$\left[\mathregular{Ca}^{2+}\right]_\mathregular{in}$ (mM)',fontsize=14)
ax5.set_xlabel(r'$t$ (ms)',fontsize=14)
ax5.set_ylabel(r'$I_\mathregular{Ca}$ (nA)',fontsize=14)
ax6.set_xlabel(r'$t$ (ms)',fontsize=14)
ax6.set_ylabel(r'$g_\mathregular{Ca}$ (mV)',fontsize=14)
ax7.set_xlabel(r'$t$ (ms)',fontsize=14)
ax7.set_ylabel(r'$I_\mathregular{SK}$ (nA)',fontsize=14)
ax8.set_xlabel(r'$t$ (ms)',fontsize=14)
ax8.set_ylabel(r'$g_\mathregular{SK}/\bar{g}_\mathregular{SK}$',fontsize=14)


N = len(t_sk1)

commonshift = -300
ishift0 = -950
ishift1 = -25+commonshift
ishift2 = -38+commonshift
ishift3 = 1999
ishift4 = 267
ishift5 = 3973

totalshift = 70
shorten1   = 100#72#989
istart     = int(9*N/16)+1550
di         = (int(10*N/16)-2850-istart)/5.
iend       = istart+int(di)
istart     = istart+shorten1+totalshift+100
iend       = iend-shorten1+totalshift+1

istart0    = istart-ishift0
istart1    = istart-ishift1
istart2    = istart-ishift2
istart3    = istart-ishift3
istart4    = istart-ishift4
istart5    = istart-ishift5
iend0      = iend-ishift0
iend1      = iend-ishift1 
iend2      = iend-ishift2 
iend3      = iend-ishift3 
iend4      = iend-ishift4 
iend5      = iend-ishift5


t1_shift_sk = t_nosk[istart0+ishift0:iend0+ishift0]
V1_shifted_CaHVA_SK = V_nosk[istart0:iend0]
V2_shifted_CaHVA_SK = V_sk1[istart1:iend1]
V3_shifted_CaHVA_SK = V_sk2[istart2:iend2]
V4_shifted_CaHVA_SK = V_sk3[istart3:iend3]
V5_shifted_CaHVA_SK = V_sk4[istart4:iend4]
V6_shifted_CaHVA_SK = V_sk5[istart5:iend5]


eca_sk1_shifted_CaHVA_SK = eca_nosk[istart0:iend0]
eca_sk2_shifted_CaHVA_SK = eca_sk1[istart1:iend1]
eca_sk3_shifted_CaHVA_SK = eca_sk2[istart2:iend2]
eca_sk4_shifted_CaHVA_SK = eca_sk3[istart3:iend3]
eca_sk5_shifted_CaHVA_SK = eca_sk4[istart4:iend4]
eca_sk6_shifted_CaHVA_SK = eca_sk5[istart5:iend5]

cai_sk1_shifted_CaHVA_SK = cai_nosk[istart0:iend0]
cai_sk2_shifted_CaHVA_SK = cai_sk1[istart1:iend1]
cai_sk3_shifted_CaHVA_SK = cai_sk2[istart2:iend2]
cai_sk4_shifted_CaHVA_SK = cai_sk3[istart3:iend3]
cai_sk5_shifted_CaHVA_SK = cai_sk4[istart4:iend4]
cai_sk6_shifted_CaHVA_SK = cai_sk5[istart5:iend5]

I_SK_sk1_shifted_CaHVA_SK = I_SK_nosk[istart0:iend0]
I_SK_sk2_shifted_CaHVA_SK = I_SK_sk1[istart1:iend1]
I_SK_sk3_shifted_CaHVA_SK = I_SK_sk2[istart2:iend2]
I_SK_sk4_shifted_CaHVA_SK = I_SK_sk3[istart3:iend3]
I_SK_sk5_shifted_CaHVA_SK = I_SK_sk4[istart4:iend4]
I_SK_sk6_shifted_CaHVA_SK = I_SK_sk5[istart5:iend5]


I_Ca_HVA_sk1_shifted_CaHVA_SK = I_Ca_HVA_nosk[istart0:iend0]
I_Ca_HVA_sk2_shifted_CaHVA_SK = I_Ca_HVA_sk1[istart1:iend1]
I_Ca_HVA_sk3_shifted_CaHVA_SK = I_Ca_HVA_sk2[istart2:iend2]
I_Ca_HVA_sk4_shifted_CaHVA_SK = I_Ca_HVA_sk3[istart3:iend3]
I_Ca_HVA_sk5_shifted_CaHVA_SK = I_Ca_HVA_sk4[istart4:iend4]
I_Ca_HVA_sk6_shifted_CaHVA_SK = I_Ca_HVA_sk5[istart5:iend5]

g_SK_sk1_shifted_CaHVA_SK = g_SK_nosk[istart0:iend0]
g_SK_sk2_shifted_CaHVA_SK = g_SK_sk1[istart1:iend1]
g_SK_sk3_shifted_CaHVA_SK = g_SK_sk2[istart2:iend2]
g_SK_sk4_shifted_CaHVA_SK = g_SK_sk3[istart3:iend3]
g_SK_sk5_shifted_CaHVA_SK = g_SK_sk4[istart4:iend4]
g_SK_sk6_shifted_CaHVA_SK = g_SK_sk5[istart5:iend5]

g_Ca_HVA_sk1_shifted_CaHVA_SK = g_Ca_HVA_nosk[istart0:iend0]
g_Ca_HVA_sk2_shifted_CaHVA_SK = g_Ca_HVA_sk1[istart1:iend1]
g_Ca_HVA_sk3_shifted_CaHVA_SK = g_Ca_HVA_sk2[istart2:iend2]
g_Ca_HVA_sk4_shifted_CaHVA_SK = g_Ca_HVA_sk3[istart3:iend3]
g_Ca_HVA_sk5_shifted_CaHVA_SK = g_Ca_HVA_sk4[istart4:iend4]
g_Ca_HVA_sk6_shifted_CaHVA_SK = g_Ca_HVA_sk5[istart5:iend5]


print('len(V1_shifted_CaHVA_SK):',len(V1_shifted_CaHVA_SK))

deltaV1 = V2_shifted_CaHVA_SK-V1_shifted_CaHVA_SK
deltaV2 = V3_shifted_CaHVA_SK-V1_shifted_CaHVA_SK
deltaV3 = V4_shifted_CaHVA_SK-V1_shifted_CaHVA_SK
deltaV4 = V5_shifted_CaHVA_SK-V1_shifted_CaHVA_SK
deltaV5 = V6_shifted_CaHVA_SK-V1_shifted_CaHVA_SK


ax1.plot(t_nosk[istart0+ishift0:iend0+ishift0], V1_shifted_CaHVA_SK,label=r'0$\bar{g}_\mathregular{SK}$')
ax1.plot(t_sk1[istart1+ishift1:iend1+ishift1], V2_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk))
ax1.plot(t_sk2[istart2+ishift2:iend2+ishift2], V3_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk2))
ax1.plot(t_sk3[istart3+ishift3:iend3+ishift3], V4_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk3))
ax1.plot(t_sk4[istart4+ishift4:iend4+ishift4], V5_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk4))
ax1.plot(t_sk5[istart5+ishift5:iend5+ishift5], V6_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk5))
ax1.set_xlim(left=t_nosk[istart0+ishift0],right=t_nosk[iend0+ishift0])
#ax2.axis([634.84,671.6,-78.2,51.1])
ax1.legend(loc='upper left',ncol=1,fontsize=11)

ax2.plot(t_nosk[istart0+ishift0:iend0+ishift0], eca_sk1_shifted_CaHVA_SK,label=r'0$\bar{g}_\mathregular{SK}$')
ax2.plot(t_sk1[istart1+ishift1:iend1+ishift1], eca_sk2_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk))
ax2.plot(t_sk2[istart2+ishift2:iend2+ishift2], eca_sk3_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk2))
ax2.plot(t_sk3[istart3+ishift3:iend3+ishift3], eca_sk4_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk3))
ax2.plot(t_sk4[istart4+ishift4:iend4+ishift4], eca_sk5_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk4))
ax2.plot(t_sk5[istart5+ishift5:iend5+ishift5], eca_sk6_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk5))
ax2.set_xlim(left=t_nosk[istart0+ishift0],right=t_nosk[iend0+ishift0])
#ax2.axis([634.84,671.6,-78.2,51.1])
ax2.legend(loc='upper right',ncol=1,fontsize=11)


ax3.plot(t1_shift_sk, deltaV1,label=r'$V_{%s\bar{g}_\mathregular{SK}}$-$V_{0\bar{g}_\mathregular{SK}}$' % str(gsk),color=mycolors[0])
ax3.plot(t1_shift_sk, deltaV2,label=r'$V_{%s\bar{g}_\mathregular{SK}}$-$V_{0\bar{g}_\mathregular{SK}}$' % str(gsk2),color=mycolors[1])
ax3.plot(t1_shift_sk, deltaV3,label=r'$V_{%s\bar{g}_\mathregular{SK}}$-$V_{0\bar{g}_\mathregular{SK}}$' % str(gsk3),color=mycolors[2])
ax3.plot(t1_shift_sk, deltaV4,label=r'$V_{%s\bar{g}_\mathregular{SK}}$-$V_{0\bar{g}_\mathregular{SK}}$' % str(gsk4),color=mycolors[3])
ax3.plot(t1_shift_sk, deltaV5,label=r'$V_{%s\bar{g}_\mathregular{SK}}$-$V_{0\bar{g}_\mathregular{SK}}$' % str(gsk5),color=mycolors[4])
ax3.set_xlim(left=t1_shift_sk[0],right=t1_shift_sk[-1])
#ax3.axis([720.109,749.359,-78.2,51.1])
ax3.legend(loc='lower left',ncol=1,fontsize=11)

ax4.plot(t_nosk[istart0+ishift0:iend0+ishift0], cai_sk1_shifted_CaHVA_SK,label=r'0$\bar{g}_\mathregular{SK}$')
ax4.plot(t_sk1[istart1+ishift1:iend1+ishift1], cai_sk2_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk))
ax4.plot(t_sk2[istart2+ishift2:iend2+ishift2], cai_sk3_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk2))
ax4.plot(t_sk3[istart3+ishift3:iend3+ishift3], cai_sk4_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk3))
ax4.plot(t_sk4[istart4+ishift4:iend4+ishift4], cai_sk5_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk4))
ax4.plot(t_sk5[istart5+ishift5:iend5+ishift5], cai_sk6_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk5))
ax4.set_xlim(left=t_nosk[istart0+ishift0],right=t_nosk[iend0+ishift0])
#ax2.axis([634.84,671.6,-78.2,51.1])
ax4.legend(loc='upper left',ncol=1,fontsize=11)


ax5.plot(t_nosk[istart0+ishift0:iend0+ishift0], I_Ca_HVA_sk1_shifted_CaHVA_SK,label=r'0$\bar{g}_\mathregular{SK}$')
ax5.plot(t_sk1[istart1+ishift1:iend1+ishift1], I_Ca_HVA_sk2_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk))
ax5.plot(t_sk2[istart2+ishift2:iend2+ishift2], I_Ca_HVA_sk3_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk2))
ax5.plot(t_sk3[istart3+ishift3:iend3+ishift3], I_Ca_HVA_sk4_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk3))
ax5.plot(t_sk4[istart4+ishift4:iend4+ishift4], I_Ca_HVA_sk5_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk4))
ax5.plot(t_sk5[istart5+ishift5:iend5+ishift5], I_Ca_HVA_sk6_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk5))
ax5.set_xlim(left=t_nosk[istart0+ishift0],right=t_nosk[iend0+ishift0])
#ax2.axis([634.84,671.6,-78.2,51.1])
ax5.legend(loc='lower left',ncol=1,fontsize=11)

ax6.plot(t_nosk[istart0+ishift0:iend0+ishift0], g_Ca_HVA_sk1_shifted_CaHVA_SK,label=r'0$\bar{g}_\mathregular{SK}$')
ax6.plot(t_sk1[istart1+ishift1:iend1+ishift1], g_Ca_HVA_sk2_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk))
ax6.plot(t_sk2[istart2+ishift2:iend2+ishift2], g_Ca_HVA_sk3_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk2))
ax6.plot(t_sk3[istart3+ishift3:iend3+ishift3], g_Ca_HVA_sk4_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk3))
ax6.plot(t_sk4[istart4+ishift4:iend4+ishift4], g_Ca_HVA_sk5_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk4))
ax6.plot(t_sk5[istart5+ishift5:iend5+ishift5], g_Ca_HVA_sk6_shifted_CaHVA_SK,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk5))
ax6.set_xlim(left=t_nosk[istart0+ishift0],right=t_nosk[iend0+ishift0])
#ax2.axis([634.84,671.6,-78.2,51.1])
ax6.legend(loc='upper left',ncol=1,fontsize=11)

#ax7.plot(t_nosk[istart:iend], I_SK_sk1_shifted_CaHVA_SK,label=r'0$\bar{g}_\mathregular{SK}$')
ax7.plot(t_sk1[istart1+ishift1:iend1+ishift1], I_SK_sk2_shifted_CaHVA_SK,color='tab:orange',label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk))
ax7.plot(t_sk2[istart2+ishift2:iend2+ishift2], I_SK_sk3_shifted_CaHVA_SK,color='tab:green',label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk2))
ax7.plot(t_sk3[istart3+ishift3:iend3+ishift3], I_SK_sk4_shifted_CaHVA_SK,color='tab:red',label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk3))
ax7.plot(t_sk4[istart4+ishift4:iend4+ishift4], I_SK_sk5_shifted_CaHVA_SK,color='tab:purple',label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk4))
ax7.plot(t_sk5[istart5+ishift5:iend5+ishift5], I_SK_sk6_shifted_CaHVA_SK,color='tab:brown',label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk5))
ax7.set_xlim(left=t_nosk[istart0+ishift0],right=t_nosk[iend0+ishift0])
#ax2.axis([634.84,671.6,-78.2,51.1])
ax7.legend(loc='upper left',ncol=1,fontsize=11)

#ax8.plot(t_nosk[istart:iend], g_SK_sk1_shifted_CaHVA_SK,label=r'0$\bar{g}_\mathregular{SK}$')
ax8.plot(t_sk1[istart1+ishift1:iend1+ishift1], g_SK_sk2_shifted_CaHVA_SK,color='tab:orange',label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk))
ax8.plot(t_sk2[istart2+ishift2:iend2+ishift2], g_SK_sk3_shifted_CaHVA_SK,color='tab:green',label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk2))
ax8.plot(t_sk3[istart3+ishift3:iend3+ishift3], g_SK_sk4_shifted_CaHVA_SK,color='tab:red',label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk3))
ax8.plot(t_sk4[istart4+ishift4:iend4+ishift4], g_SK_sk5_shifted_CaHVA_SK,color='tab:purple',label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk4))
ax8.plot(t_sk5[istart5+ishift5:iend5+ishift5], g_SK_sk6_shifted_CaHVA_SK,color='tab:brown',label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk5))
ax8.set_xlim(left=t_nosk[istart0+ishift0],right=t_nosk[iend0+ishift0])
#ax2.axis([634.84,671.6,-78.2,51.1])
ax8.legend(loc='upper center',ncol=1,fontsize=11)

'''
print('len(t1_shift):',len(t1_shift))
print('len(t1_shift_sk):',len(t1_shift_sk))
'''

fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig('Results/Soma10/Compare/trace_cainfo_idur'+str(idur)+'_iamp'+str(iamp)+'_SK_aligned.png')
plt.show()

