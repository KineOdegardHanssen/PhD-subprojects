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
	
startconcentration = 0.0001
halfconcentration   = 0.00043287612810830614

idur = 1000
iamp = 0.02

plotitall = False # True # 

gcahva = 0.2
gcahva_label = 1.0
gcahva_base  = 0.00014931667074610222
Vshifts = [-20,-10,0,10,20]
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

t1 = ts[0] # 0.2
t2 = ts[1] # 0.5
t3 = ts[2] # 1.0
t4 = ts[3] # 1.5
t5 = ts[4] # 2.0
V1 = Vs[0]
V2 = Vs[1]
V3 = Vs[2]
V4 = Vs[3]
V5 = Vs[4]
eca1 = ecas[0]
eca2 = ecas[1]
eca3 = ecas[2]
eca4 = ecas[3]
eca5 = ecas[4]
###
cai1 = cais[0]
cai2 = cais[1]
cai3 = cais[2]
cai4 = cais[3]
cai5 = cais[4]
###
cao1 = caos[0]
cao2 = caos[1]
cao3 = caos[2]
cao4 = caos[3]
cao5 = caos[4]
###
I_Ca_HVA1 = I_Ca_HVAs[0]
I_Ca_HVA2 = I_Ca_HVAs[1]
I_Ca_HVA3 = I_Ca_HVAs[2]
I_Ca_HVA4 = I_Ca_HVAs[3]
I_Ca_HVA5 = I_Ca_HVAs[4]
###
g_Ca_HVA1 = g_Ca_HVAs[0]
g_Ca_HVA2 = g_Ca_HVAs[1]
g_Ca_HVA3 = g_Ca_HVAs[2]
g_Ca_HVA4 = g_Ca_HVAs[3]
g_Ca_HVA5 = g_Ca_HVAs[4]



ishift1 = 0#-200#2426
ishift2 = 921#1069
ishift3 = -180
ishift4 = 5422
ishift5 = 995

shiftall = -100#200
spdelta = 2200
epdelta = 2200
# Don't know these, should test out:
istart = int(9*N/16)+3938+shiftall+spdelta+100
iend   = int(5*N/8)-3848+shiftall+epdelta
istart1 = istart-ishift1
istart2 = istart-ishift2
istart3 = istart-ishift3
istart4 = istart-ishift4
istart5 = istart-ishift5
iend1   = iend-ishift1
iend2   = iend-ishift2
iend3   = iend-ishift3
iend4   = iend-ishift4
iend5   = iend-ishift5

t1_shift = t1[istart1+ishift1:iend1+ishift1]
t2_shift = t2[istart2+ishift2:iend2+ishift2]
t3_shift = t3[istart3+ishift3:iend3+ishift3]
t4_shift = t4[istart4+ishift4:iend4+ishift4]
t5_shift = t5[istart5+ishift5:iend5+ishift5]

V1_shifted_CaHVA = V1[istart1:iend1]
V2_shifted_CaHVA = V2[istart2:iend2]
V3_shifted_CaHVA = V3[istart3:iend3]
V4_shifted_CaHVA = V4[istart4:iend4]
V5_shifted_CaHVA = V5[istart5:iend5]

eca1_shifted_CaHVA = eca1[istart1:iend1]
eca2_shifted_CaHVA = eca2[istart2:iend2]
eca3_shifted_CaHVA = eca3[istart3:iend3] 
eca4_shifted_CaHVA = eca4[istart4:iend4] 
eca5_shifted_CaHVA = eca5[istart5:iend5]
###
cai1_shifted_CaHVA = cai1[istart1:iend1]
cai2_shifted_CaHVA = cai2[istart2:iend2]
cai3_shifted_CaHVA = cai3[istart3:iend3]
cai4_shifted_CaHVA = cai4[istart4:iend4]
cai5_shifted_CaHVA = cai5[istart5:iend5]
###
I_Ca_HVA1_shifted_CaHVA = I_Ca_HVA1[istart1:iend1]
I_Ca_HVA2_shifted_CaHVA = I_Ca_HVA2[istart2:iend2]
I_Ca_HVA3_shifted_CaHVA = I_Ca_HVA3[istart3:iend3]
I_Ca_HVA4_shifted_CaHVA = I_Ca_HVA4[istart4:iend4]
I_Ca_HVA5_shifted_CaHVA = I_Ca_HVA5[istart5:iend5]
###
g_Ca_HVA1_shifted_CaHVA = g_Ca_HVA1[istart1:iend1]
g_Ca_HVA2_shifted_CaHVA = g_Ca_HVA2[istart2:iend2]
g_Ca_HVA3_shifted_CaHVA = g_Ca_HVA3[istart3:iend3]
g_Ca_HVA4_shifted_CaHVA = g_Ca_HVA4[istart4:iend4]
g_Ca_HVA5_shifted_CaHVA = g_Ca_HVA5[istart5:iend5]

deltaV1 = V1_shifted_CaHVA-V3_shifted_CaHVA
deltaV2 = V2_shifted_CaHVA-V3_shifted_CaHVA
deltaV3 = V4_shifted_CaHVA-V3_shifted_CaHVA
deltaV4 = V5_shifted_CaHVA-V3_shifted_CaHVA

print('len(V1_shifted_CaHVA):',len(V1_shifted_CaHVA))
print('len(V2_shifted_CaHVA):',len(V2_shifted_CaHVA))
print('len(V3_shifted_CaHVA):',len(V3_shifted_CaHVA))
print('len(V4_shifted_CaHVA):',len(V4_shifted_CaHVA))
#print('len(V5_shifted_CaHVA):',len(V5_shifted_CaHVA))


peaktimes1 = givepeak(t1_shift,V1_shifted_CaHVA)
peaktimes2 = givepeak(t2_shift,V2_shifted_CaHVA)
peaktimes3 = givepeak(t3_shift,V3_shifted_CaHVA)
peaktimes4 = givepeak(t4_shift,V4_shifted_CaHVA)
peaktimes5 = givepeak(t5_shift,V5_shifted_CaHVA)


peaktime_CaHVA = peaktimes1[0]

if plotitall==True:
    plt.figure()
    plt.plot(t5[istart5+ishift5:iend5+ishift5], V5_shifted_CaHVA,label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
    plt.plot(t4[istart4+ishift4:iend4+ishift4], V4_shifted_CaHVA,label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
    plt.plot(t3[istart3+ishift3:iend3+ishift3], V3_shifted_CaHVA,label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
    plt.plot(t2[istart2+ishift2:iend2+ishift2], V2_shifted_CaHVA,label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
    plt.plot(t1[istart1+ishift1:iend1+ishift1], V1_shifted_CaHVA,label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
    plt.legend(loc='upper left')
    plt.xlabel(r'$t$ (s)')
    plt.ylabel(r'$V$ (mV)')
    plt.show()


plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

fig = plt.figure(figsize=(10,10),dpi=300)#(figsize=(8,3),dpi=300)

gs = gridspec.GridSpec(3, 4)
    
ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])    
ax3 = plt.subplot(gs[1, 0:2])
ax4 = plt.subplot(gs[1, 2:4])    
ax5 = plt.subplot(gs[2, 0:2])
ax6 = plt.subplot(gs[2, 2:4])
    
#fig.suptitle(r'$I=$%s nA' % str(iamp),fontsize=20)

ax1.set_title(r'A',loc='left',fontsize=14)
ax2.set_title(r'B',loc='left',fontsize=14)
ax3.set_title(r'C',loc='left',fontsize=14)
ax4.set_title(r'D',loc='left',fontsize=14)
ax5.set_title(r'E',loc='left',fontsize=14)
ax6.set_title(r'F',loc='left',fontsize=14)
ax1.set_title(r'Voltage trace',fontsize=14)#loc='right',
ax2.set_title(r'$E_\mathregular{Ca}$',fontsize=14)#loc='right',
ax3.set_title(r'Voltage trace difference',fontsize=14)#loc='right',
ax4.set_title(r'$\left[\mathregular{Ca}^{2+}\right]_\mathregular{in}$',fontsize=14)#loc='right',
ax5.set_title(r'$I_\mathregular{Ca}$',fontsize=14)#loc='right',
ax6.set_title(r'$g_\mathregular{Ca}/\bar{g}_\mathregular{Ca}$',fontsize=14)#loc='right',
mycolors = ['tab:orange','tab:green','tab:red','tab:purple','tab:brown']
    
ax1.set_xlabel(r'$t$ (ms)',fontsize=12)
ax1.set_ylabel(r'$V$ (mV)',fontsize=12)
ax2.set_xlabel(r'$t$ (ms)',fontsize=12)
ax2.set_ylabel(r'$V$ (mV)',fontsize=12)
ax3.set_xlabel(r'$t$ (ms)',fontsize=12)
ax3.set_ylabel(r'$\Delta V$ (mV)',fontsize=12)
ax4.set_xlabel(r'$t$ (ms)',fontsize=12)
ax4.set_ylabel(r'$\left[\mathregular{Ca}^{2+}\right]_\mathregular{in}$ (mM)',fontsize=12)
ax5.set_xlabel(r'$t$ (ms)',fontsize=12)
ax5.set_ylabel(r'$I_\mathregular{Ca}$ (nA)',fontsize=12)
ax6.set_xlabel(r'$t$ (ms)',fontsize=12)
ax6.set_ylabel(r'$g_\mathregular{Ca}/\bar{g}_\mathregular{Ca}$',fontsize=12)

#'''

ax1.axvline(x=peaktimes1,color='k',linestyle='--',linewidth=0.75)
ax1.plot(t5[istart5+ishift5:iend5+ishift5], V5_shifted_CaHVA,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax1.plot(t4[istart4+ishift4:iend4+ishift4], V4_shifted_CaHVA,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax1.plot(t3[istart3+ishift3:iend3+ishift3], V3_shifted_CaHVA,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax1.plot(t2[istart2+ishift2:iend2+ishift2], V2_shifted_CaHVA,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax1.plot(t1[istart1+ishift1:iend1+ishift1], V1_shifted_CaHVA,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax1.set_xlim(left=t1[istart1+ishift1],right=t1[iend1+ishift1])
ax1.legend(loc='upper right',ncol=1,fontsize=10)

ax2.axvline(x=peaktimes1,color='k',linestyle='--',linewidth=0.75)
ax2.plot(t5[istart5+ishift5:iend5+ishift5], eca5_shifted_CaHVA,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax2.plot(t4[istart4+ishift4:iend4+ishift4], eca4_shifted_CaHVA,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax2.plot(t3[istart3+ishift3:iend3+ishift3], eca3_shifted_CaHVA,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax2.plot(t2[istart2+ishift2:iend2+ishift2], eca2_shifted_CaHVA,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax2.plot(t1[istart1+ishift1:iend1+ishift1], eca1_shifted_CaHVA,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax2.set_xlim(left=t1[istart1+ishift1],right=t1[iend1+ishift1])
#ax2.legend(loc='upper left',ncol=1,fontsize=11)

ax3.axvline(x=peaktimes1,color='k',linestyle='--',linewidth=0.75)
ax3.plot(t1_shift, deltaV4,color=colors[0],label=r'$V_{%i\ \mathregular{mV}}-V_{0\ \mathregular{mV}}$' % Vshifts[0])
ax3.plot(t1_shift, deltaV3,color=colors[1],label=r'$V_{%i\ \mathregular{mV}}-V_{0\ \mathregular{mV}}$' % Vshifts[1])
ax3.plot(t1_shift, deltaV2,color=colors[3],label=r'$V_{%i\ \mathregular{mV}}-V_{0\ \mathregular{mV}}$' % Vshifts[3])
ax3.plot(t1_shift, deltaV1,color=colors[4],label=r'$V_{%i\ \mathregular{mV}}-V_{0\ \mathregular{mV}}$' % Vshifts[4])
ax3.set_xlim(left=t1[istart1+ishift1],right=t1[iend1+ishift1])
ax3.legend(loc='upper right',ncol=1,fontsize=11)

ax4.axvline(x=peaktimes1,color='k',linestyle='--',linewidth=0.75)
ax4.axhline(y=startconcentration,color='grey',linestyle=':',linewidth=0.75)
ax4.axhline(y=halfconcentration,color='k',linestyle='--',linewidth=0.75)
ax4.plot(t5[istart5+ishift5:iend5+ishift5], cai5_shifted_CaHVA,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax4.plot(t4[istart4+ishift4:iend4+ishift4], cai4_shifted_CaHVA,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax4.plot(t3[istart3+ishift3:iend3+ishift3], cai3_shifted_CaHVA,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax4.plot(t2[istart2+ishift2:iend2+ishift2], cai2_shifted_CaHVA,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax4.plot(t1[istart1+ishift1:iend1+ishift1], cai1_shifted_CaHVA,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax4.set_xlim(left=t1[istart1+ishift1],right=t1[iend1+ishift1])
#ax4.legend(loc='upper left',ncol=1,fontsize=11)

ax5.axvline(x=peaktimes1,color='k',linestyle='--',linewidth=0.75)
ax5.plot(t5[istart5+ishift5:iend5+ishift5], I_Ca_HVA5_shifted_CaHVA,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax5.plot(t4[istart4+ishift4:iend4+ishift4], I_Ca_HVA4_shifted_CaHVA,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax5.plot(t3[istart3+ishift3:iend3+ishift3], I_Ca_HVA3_shifted_CaHVA,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax5.plot(t2[istart2+ishift2:iend2+ishift2], I_Ca_HVA2_shifted_CaHVA,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax5.plot(t1[istart1+ishift1:iend1+ishift1], I_Ca_HVA1_shifted_CaHVA,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax5.set_xlim(left=t1[istart1+ishift1],right=t1[iend1+ishift1])
#ax5.legend(loc='lower left',ncol=1,fontsize=10)

ax6.axvline(x=peaktimes1,color='k',linestyle='--',linewidth=0.75)
ax6.plot(t5[istart5+ishift5:iend5+ishift5], g_Ca_HVA5_shifted_CaHVA,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax6.plot(t4[istart4+ishift4:iend4+ishift4], g_Ca_HVA4_shifted_CaHVA,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax6.plot(t3[istart3+ishift3:iend3+ishift3], g_Ca_HVA3_shifted_CaHVA,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax6.plot(t2[istart2+ishift2:iend2+ishift2], g_Ca_HVA2_shifted_CaHVA,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax6.plot(t1[istart1+ishift1:iend1+ishift1], g_Ca_HVA1_shifted_CaHVA,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax6.set_xlim(left=t1[istart1+ishift1],right=t1[iend1+ishift1])
#ax6.legend(loc='upper left',ncol=1,fontsize=10)

fig.tight_layout()#(rect=[0, 0.03, 1, 0.95])
plt.savefig('Results/Soma10/Compare/trace_cainfo_idur'+str(idur)+'_iamp'+str(iamp)+'_CaHVA_SHIFTS_aligned.png')
#'''

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

t_sk1 = ts[0] # 0.2
t_sk2 = ts[1] # 0.5
t_sk3 = ts[2] # 1.0
t_sk4 = ts[3] # 1.5
t_sk5 = ts[4] # 2.0
V_sk1 = Vs[0]
V_sk2 = Vs[1]
V_sk3 = Vs[2]
V_sk4 = Vs[3]
V_sk5 = Vs[4]


eca_sk1 = ecas[0]
eca_sk2 = ecas[1]
eca_sk3 = ecas[2]
eca_sk4 = ecas[3]
eca_sk5 = ecas[4]

cai_sk1 = cais[0]
cai_sk2 = cais[1]
cai_sk3 = cais[2]
cai_sk4 = cais[3]
cai_sk5 = cais[4]

cao_sk1 = caos[0]
cao_sk2 = caos[1]
cao_sk3 = caos[2]
cao_sk4 = caos[3]
cao_sk5 = caos[4]

I_SK_sk1 = I_SKs[0]
I_SK_sk2 = I_SKs[1]
I_SK_sk3 = I_SKs[2]
I_SK_sk4 = I_SKs[3]
I_SK_sk5 = I_SKs[4]

g_SK_sk1 = g_SKs[0]
g_SK_sk2 = g_SKs[1]
g_SK_sk3 = g_SKs[2]
g_SK_sk4 = g_SKs[3]
g_SK_sk5 = g_SKs[4]

I_Ca_HVA_sk1 = I_Ca_HVAs[0]
I_Ca_HVA_sk2 = I_Ca_HVAs[1]
I_Ca_HVA_sk3 = I_Ca_HVAs[2]
I_Ca_HVA_sk4 = I_Ca_HVAs[3]
I_Ca_HVA_sk5 = I_Ca_HVAs[4]

g_Ca_HVA_sk1 = g_Ca_HVAs[0]
g_Ca_HVA_sk2 = g_Ca_HVAs[1]
g_Ca_HVA_sk3 = g_Ca_HVAs[2]
g_Ca_HVA_sk4 = g_Ca_HVAs[3]
g_Ca_HVA_sk5 = g_Ca_HVAs[4]

gsk_base = 0.0028175455472127641
g_SK_sk1 = g_SK_sk1/float(gsk*gsk_base)
g_SK_sk2 = g_SK_sk2/float(gsk*gsk_base)
g_SK_sk3 = g_SK_sk3/float(gsk*gsk_base)
g_SK_sk4 = g_SK_sk4/float(gsk*gsk_base)
g_SK_sk5 = g_SK_sk5/float(gsk*gsk_base)


g_Ca_HVA_sk1 = g_Ca_HVA_sk1/float(gcahva*gcahva_base)
g_Ca_HVA_sk2 = g_Ca_HVA_sk2/float(gcahva*gcahva_base)
g_Ca_HVA_sk3 = g_Ca_HVA_sk3/float(gcahva*gcahva_base)
g_Ca_HVA_sk4 = g_Ca_HVA_sk4/float(gcahva*gcahva_base)
g_Ca_HVA_sk5 = g_Ca_HVA_sk5/float(gcahva*gcahva_base)

Nt = len(t_sk1)


ishift1 = 150 #-300
ishift2 = -225 #-213
ishift3 = 142 #-242
ishift4 = -717 #-194
ishift5 = 63 #-1262

shorten1 = 0#467
istart = int(9*N/16)+800
iend   = int(5*N/8)-7100
istart  = istart+shorten1+100
iend    = iend-shorten1-1

istart1 = istart-ishift1
istart2 = istart-ishift2
istart3 = istart-ishift3
istart4 = istart-ishift4
istart5 = istart-ishift5
iend1   = iend-ishift1
iend2   = iend-ishift2
iend3   = iend-ishift3
iend4   = iend-ishift4
iend5   = iend-ishift5


print('t_sk1[istart]:',t_sk1[istart])
print('t_sk1[iend]:',t_sk1[iend])
print('t_sk1[istart2]:',t_sk1[istart2])
print('t_sk1[iend2]:',t_sk1[iend2])
print('V_sk3[istart2]:',V_sk3[istart2])
print('V_sk3[iend2]:',V_sk3[iend2])

t1_shift_sk = t_sk1[istart1+ishift1:iend1+ishift1]
t2_shift_sk = t_sk2[istart2+ishift2:iend2+ishift2]
t3_shift_sk = t_sk3[istart3+ishift3:iend3+ishift3]
t4_shift_sk = t_sk4[istart4+ishift4:iend4+ishift4]
t5_shift_sk = t_sk5[istart5+ishift5:iend5+ishift5]
V1_shifted_CaHVA_SK = V_sk1[istart1:iend1]
V2_shifted_CaHVA_SK = V_sk2[istart2:iend2]
V3_shifted_CaHVA_SK = V_sk3[istart3:iend3]
V4_shifted_CaHVA_SK = V_sk4[istart4:iend4]
V5_shifted_CaHVA_SK = V_sk5[istart5:iend5]


peakinds1, peakvals1 = find_peakind(V1_shifted_CaHVA_SK)
peakinds2, peakvals2 = find_peakind(V2_shifted_CaHVA_SK)
peakinds3, peakvals3 = find_peakind(V3_shifted_CaHVA_SK)
peakinds4, peakvals4 = find_peakind(V4_shifted_CaHVA_SK)
peakinds5, peakvals5 = find_peakind(V5_shifted_CaHVA_SK)

peaktimes1_SK = givepeak(t1_shift_sk,V1_shifted_CaHVA_SK)
peaktimes2_SK = givepeak(t2_shift_sk,V2_shifted_CaHVA_SK)
peaktimes3_SK = givepeak(t3_shift_sk,V3_shifted_CaHVA_SK)
peaktimes4_SK = givepeak(t4_shift_sk,V4_shifted_CaHVA_SK)
peaktimes5_SK = givepeak(t5_shift_sk,V5_shifted_CaHVA_SK)

print('peaktimes1:',peaktimes1)
print('peaktimes2:',peaktimes2)
print('peaktimes3:',peaktimes3)
print('peaktimes4:',peaktimes4)
print('peaktimes5:',peaktimes5)

peaktime_CaHVA_SK = peaktimes1_SK[0]

if plotitall==True:
    plt.figure()
    plt.plot(t_sk5[istart+ishift5:iend+ishift5], V5_shifted_CaHVA_SK,label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
    plt.plot(t_sk4[istart+ishift4:iend+ishift4], V4_shifted_CaHVA_SK,label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
    plt.plot(t_sk3[istart+ishift3:iend+ishift3], V3_shifted_CaHVA_SK,label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
    plt.plot(t_sk2[istart+ishift2:iend+ishift2], V2_shifted_CaHVA_SK,label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
    plt.plot(t_sk1[istart+ishift1:iend+ishift1], V1_shifted_CaHVA_SK,label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
    plt.legend(loc='upper left')
    plt.xlabel(r'$t$ (s)')
    plt.ylabel(r'$V$ (mV)')
    plt.show()
'''
print('peakinds1:', peakinds1)
print('peakinds2:', peakinds2)
print('peakinds3:', peakinds3)
print('peakinds4:', peakinds4)
print('peakinds5:', peakinds5)
'''

eca_sk1_shifted_CaHVA_SK = eca_sk1[istart1:iend1]
eca_sk2_shifted_CaHVA_SK = eca_sk2[istart2:iend2]
eca_sk3_shifted_CaHVA_SK = eca_sk3[istart3:iend3]
eca_sk4_shifted_CaHVA_SK = eca_sk4[istart4:iend4]
eca_sk5_shifted_CaHVA_SK = eca_sk5[istart5:iend5]

cai_sk1_shifted_CaHVA_SK = cai_sk1[istart1:iend1]
cai_sk2_shifted_CaHVA_SK = cai_sk2[istart2:iend2]
cai_sk3_shifted_CaHVA_SK = cai_sk3[istart3:iend3]
cai_sk4_shifted_CaHVA_SK = cai_sk4[istart4:iend4]
cai_sk5_shifted_CaHVA_SK = cai_sk5[istart5:iend5]

I_SK_sk1_shifted_CaHVA_SK = I_SK_sk1[istart1:iend1]
I_SK_sk2_shifted_CaHVA_SK = I_SK_sk2[istart2:iend2]
I_SK_sk3_shifted_CaHVA_SK = I_SK_sk3[istart3:iend3]
I_SK_sk4_shifted_CaHVA_SK = I_SK_sk4[istart4:iend4]
I_SK_sk5_shifted_CaHVA_SK = I_SK_sk5[istart5:iend5]

I_Ca_HVA_sk1_shifted_CaHVA_SK = I_Ca_HVA_sk1[istart1:iend1]
I_Ca_HVA_sk2_shifted_CaHVA_SK = I_Ca_HVA_sk2[istart2:iend2]
I_Ca_HVA_sk3_shifted_CaHVA_SK = I_Ca_HVA_sk3[istart3:iend3]
I_Ca_HVA_sk4_shifted_CaHVA_SK = I_Ca_HVA_sk4[istart4:iend4]
I_Ca_HVA_sk5_shifted_CaHVA_SK = I_Ca_HVA_sk5[istart5:iend5]

g_SK_sk1_shifted_CaHVA_SK = g_SK_sk1[istart1:iend1]
g_SK_sk2_shifted_CaHVA_SK = g_SK_sk2[istart2:iend2]
g_SK_sk3_shifted_CaHVA_SK = g_SK_sk3[istart3:iend3]
g_SK_sk4_shifted_CaHVA_SK = g_SK_sk4[istart4:iend4]
g_SK_sk5_shifted_CaHVA_SK = g_SK_sk5[istart5:iend5]

g_Ca_HVA_sk1_shifted_CaHVA_SK = g_Ca_HVA_sk1[istart1:iend1]
g_Ca_HVA_sk2_shifted_CaHVA_SK = g_Ca_HVA_sk2[istart2:iend2]
g_Ca_HVA_sk3_shifted_CaHVA_SK = g_Ca_HVA_sk3[istart3:iend3]
g_Ca_HVA_sk4_shifted_CaHVA_SK = g_Ca_HVA_sk4[istart4:iend4]
g_Ca_HVA_sk5_shifted_CaHVA_SK = g_Ca_HVA_sk5[istart5:iend5]


print('len(V1_shifted_CaHVA_SK):',len(V1_shifted_CaHVA_SK))

deltaV1_SK = V1_shifted_CaHVA_SK-V3_shifted_CaHVA_SK
deltaV2_SK = V2_shifted_CaHVA_SK-V3_shifted_CaHVA_SK
deltaV3_SK = V4_shifted_CaHVA_SK-V3_shifted_CaHVA_SK
deltaV4_SK = V5_shifted_CaHVA_SK-V3_shifted_CaHVA_SK


fig = plt.figure(figsize=(10,12),dpi=300)#(figsize=(8,3),dpi=300)

gs = gridspec.GridSpec(4, 4)
    
ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])    
ax3 = plt.subplot(gs[1, 0:2])
ax4 = plt.subplot(gs[1, 2:4])    
ax5 = plt.subplot(gs[2, 0:2])
ax6 = plt.subplot(gs[2, 2:4])    
ax7 = plt.subplot(gs[3, 0:2])
ax8 = plt.subplot(gs[3, 2:4])
    
#fig.suptitle(r'$I=$%s nA' % str(iamp),fontsize=20)

ax1.set_title(r'A',loc='left',fontsize=14)
ax2.set_title(r'B',loc='left',fontsize=14)
ax3.set_title(r'C',loc='left',fontsize=14)
ax4.set_title(r'D',loc='left',fontsize=14)
ax5.set_title(r'E',loc='left',fontsize=14)
ax6.set_title(r'F',loc='left',fontsize=14)
ax7.set_title(r'G',loc='left',fontsize=14)
ax8.set_title(r'H',loc='left',fontsize=14)
ax1.set_title(r'Voltage trace',fontsize=14)#loc='right',
ax2.set_title(r'$E_\mathregular{Ca}$',fontsize=14)#loc='right',
ax3.set_title(r'Voltage trace difference',fontsize=14)#loc='right',
ax4.set_title(r'$\left[\mathregular{Ca}^{2+}\right]_\mathregular{in}$',fontsize=14)#loc='right',
ax5.set_title(r'$I_\mathregular{Ca}$',fontsize=14)#loc='right',
ax6.set_title(r'$g_\mathregular{Ca}/\bar{g}_\mathregular{Ca}$',fontsize=14)#loc='right',
ax7.set_title(r'$I_\mathregular{SK}$',fontsize=14)#loc='right',
ax8.set_title(r'$g_\mathregular{SK}/\bar{g}_\mathregular{SK}$',fontsize=14)#loc='right',
mycolors = ['tab:orange','tab:green','tab:red','tab:purple','tab:brown']
    
ax1.set_xlabel(r'$t$ (ms)',fontsize=12)
ax1.set_ylabel(r'$V$ (mV)',fontsize=12)
ax2.set_xlabel(r'$t$ (ms)',fontsize=12)
ax2.set_ylabel(r'$V$ (mV)',fontsize=12)
ax3.set_xlabel(r'$t$ (ms)',fontsize=12)
ax3.set_ylabel(r'$\Delta V$ (mV)',fontsize=12)
ax4.set_xlabel(r'$t$ (ms)',fontsize=12)
ax4.set_ylabel(r'$\left[\mathregular{Ca}^{2+}\right]_\mathregular{in}$ (mM)',fontsize=12)
ax5.set_xlabel(r'$t$ (ms)',fontsize=12)
ax5.set_ylabel(r'$I_\mathregular{Ca}$ (nA)',fontsize=12)
ax6.set_xlabel(r'$t$ (ms)',fontsize=12)
ax6.set_ylabel(r'$g_\mathregular{Ca}/\bar{g}_\mathregular{Ca}$',fontsize=12)
ax7.set_xlabel(r'$t$ (ms)',fontsize=12)
ax7.set_ylabel(r'$I_\mathregular{SK}$ (nA)',fontsize=12)
ax8.set_xlabel(r'$t$ (ms)',fontsize=12)
ax8.set_ylabel(r'$g_\mathregular{SK}/\bar{g}_\mathregular{SK}$',fontsize=12)


ax1.axvline(x=peaktimes1_SK,color='k',linestyle='--',linewidth=0.75)
ax1.plot(t_sk5[istart5+ishift5:iend5+ishift5], V5_shifted_CaHVA_SK,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax1.plot(t_sk4[istart4+ishift4:iend4+ishift4], V4_shifted_CaHVA_SK,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax1.plot(t_sk3[istart3+ishift3:iend3+ishift3], V3_shifted_CaHVA_SK,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax1.plot(t_sk2[istart2+ishift2:iend2+ishift2], V2_shifted_CaHVA_SK,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax1.plot(t_sk1[istart1+ishift1:iend1+ishift1], V1_shifted_CaHVA_SK,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax1.set_xlim(left=t_sk1[istart1+ishift1],right=t_sk1[iend1+ishift1])
ax1.legend(loc='upper right',ncol=1,fontsize=10)

ax2.axvline(x=peaktimes1_SK,color='k',linestyle='--',linewidth=0.75)
ax2.plot(t_sk5[istart5+ishift5:iend5+ishift5], eca_sk5_shifted_CaHVA_SK,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax2.plot(t_sk4[istart4+ishift4:iend4+ishift4], eca_sk4_shifted_CaHVA_SK,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax2.plot(t_sk3[istart3+ishift3:iend3+ishift3], eca_sk3_shifted_CaHVA_SK,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax2.plot(t_sk2[istart2+ishift2:iend2+ishift2], eca_sk2_shifted_CaHVA_SK,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax2.plot(t_sk1[istart1+ishift1:iend1+ishift1], eca_sk1_shifted_CaHVA_SK,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax2.set_xlim(left=t_sk1[istart1+ishift1],right=t_sk1[iend1+ishift1])
ax2.legend(loc='lower left',ncol=1,fontsize=11)

ax3.axvline(x=peaktimes1_SK,color='k',linestyle='--',linewidth=0.75)
ax3.plot(t1_shift_sk, deltaV4_SK,color=colors[0],label=r'$V_{%i\ \mathregular{mV}}-V_{0\ \mathregular{mV}}$' % Vshifts[0])
ax3.plot(t1_shift_sk, deltaV3_SK,color=colors[1],label=r'$V_{%i\ \mathregular{mV}}-V_{0\ \mathregular{mV}}$' % Vshifts[1])
ax3.plot(t1_shift_sk, deltaV2_SK,color=colors[3],label=r'$V_{%i\ \mathregular{mV}}-V_{0\ \mathregular{mV}}$' % Vshifts[3])
ax3.plot(t1_shift_sk, deltaV1_SK,color=colors[4],label=r'$V_{%i\ \mathregular{mV}}-V_{0\ \mathregular{mV}}$' % Vshifts[4])
ax3.set_xlim(left=t1_shift_sk[0],right=t1_shift_sk[-1])
ax3.legend(loc='lower right',ncol=1,fontsize=11)

ax4.axvline(x=peaktimes1_SK,color='k',linestyle='--',linewidth=0.75)
ax4.axhline(y=startconcentration,color='grey',linestyle=':',linewidth=0.75)
ax4.axhline(y=halfconcentration,color='k',linestyle='--',linewidth=0.75)
ax4.plot(t_sk5[istart5+ishift5:iend5+ishift5], cai_sk5_shifted_CaHVA_SK,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax4.plot(t_sk4[istart4+ishift4:iend4+ishift4], cai_sk4_shifted_CaHVA_SK,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax4.plot(t_sk3[istart3+ishift3:iend3+ishift3], cai_sk3_shifted_CaHVA_SK,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax4.plot(t_sk2[istart2+ishift2:iend2+ishift2], cai_sk2_shifted_CaHVA_SK,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax4.plot(t_sk1[istart1+ishift1:iend1+ishift1], cai_sk1_shifted_CaHVA_SK,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax4.set_xlim(left=t_sk1[istart1+ishift1],right=t_sk1[iend1+ishift1])
ax4.legend(loc='upper left',ncol=1,fontsize=11)

ax5.axvline(x=peaktimes1_SK,color='k',linestyle='--',linewidth=0.75)
ax5.plot(t_sk5[istart5+ishift5:iend5+ishift5], I_Ca_HVA_sk5_shifted_CaHVA_SK,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax5.plot(t_sk4[istart4+ishift4:iend4+ishift4], I_Ca_HVA_sk4_shifted_CaHVA_SK,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax5.plot(t_sk3[istart3+ishift3:iend3+ishift3], I_Ca_HVA_sk3_shifted_CaHVA_SK,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax5.plot(t_sk2[istart2+ishift2:iend2+ishift2], I_Ca_HVA_sk2_shifted_CaHVA_SK,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax5.plot(t_sk1[istart1+ishift1:iend1+ishift1], I_Ca_HVA_sk1_shifted_CaHVA_SK,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax5.set_xlim(left=t_sk1[istart1+ishift1],right=t_sk1[iend1+ishift1])
ax5.legend(loc='lower left',ncol=1,fontsize=10)

ax6.axvline(x=peaktimes1_SK,color='k',linestyle='--',linewidth=0.75)
ax6.plot(t_sk5[istart5+ishift5:iend5+ishift5], g_Ca_HVA_sk5_shifted_CaHVA_SK,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax6.plot(t_sk4[istart4+ishift4:iend4+ishift4], g_Ca_HVA_sk4_shifted_CaHVA_SK,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax6.plot(t_sk3[istart3+ishift3:iend3+ishift3], g_Ca_HVA_sk3_shifted_CaHVA_SK,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax6.plot(t_sk2[istart2+ishift2:iend2+ishift2], g_Ca_HVA_sk2_shifted_CaHVA_SK,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax6.plot(t_sk1[istart1+ishift1:iend1+ishift1], g_Ca_HVA_sk1_shifted_CaHVA_SK,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax6.set_xlim(left=t_sk1[istart1+ishift1],right=t_sk1[iend1+ishift1])
ax6.legend(loc='upper left',ncol=1,fontsize=10)

ax7.axvline(x=peaktimes1_SK,color='k',linestyle='--',linewidth=0.75)
ax7.plot(t_sk5[istart5+ishift5:iend5+ishift5], I_SK_sk5_shifted_CaHVA_SK,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax7.plot(t_sk4[istart4+ishift4:iend4+ishift4], I_SK_sk4_shifted_CaHVA_SK,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax7.plot(t_sk3[istart3+ishift3:iend3+ishift3], I_SK_sk3_shifted_CaHVA_SK,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax7.plot(t_sk2[istart2+ishift2:iend2+ishift2], I_SK_sk2_shifted_CaHVA_SK,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax7.plot(t_sk1[istart1+ishift1:iend1+ishift1], I_SK_sk1_shifted_CaHVA_SK,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax7.set_xlim(left=t_sk1[istart1+ishift1],right=t_sk1[iend1+ishift1])
ax7.legend(loc='upper left',ncol=1,fontsize=10)

ax8.axvline(x=peaktimes1_SK,color='k',linestyle='--',linewidth=0.75)
ax8.plot(t_sk5[istart5+ishift5:iend5+ishift5], g_SK_sk5_shifted_CaHVA_SK,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax8.plot(t_sk4[istart4+ishift4:iend4+ishift4], g_SK_sk4_shifted_CaHVA_SK,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax8.plot(t_sk3[istart3+ishift3:iend3+ishift3], g_SK_sk3_shifted_CaHVA_SK,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax8.plot(t_sk2[istart2+ishift2:iend2+ishift2], g_SK_sk2_shifted_CaHVA_SK,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax8.plot(t_sk1[istart1+ishift1:iend1+ishift1], g_SK_sk1_shifted_CaHVA_SK,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax8.set_xlim(left=t_sk1[istart1+ishift1],right=t_sk1[iend1+ishift1])
ax8.legend(loc='upper left',ncol=1,fontsize=11)

'''
print('len(t1_shift):',len(t1_shift))
print('len(t1_shift_sk):',len(t1_shift_sk))
'''

fig.tight_layout()#(rect=[0, 0.03, 1, 0.95])
plt.savefig('Results/Soma10/Compare/trace_cainfo_idur'+str(idur)+'_iamp'+str(iamp)+'_SK_SHIFTS_aligned.png')
#plt.show()


plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

fig = plt.figure(figsize=(10,12),dpi=300)#(figsize=(8,3),dpi=300)

gs = gridspec.GridSpec(4, 4)
    
ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])    
ax3 = plt.subplot(gs[1, 0:2])
ax4 = plt.subplot(gs[1, 2:4])    
ax5 = plt.subplot(gs[2, 0:2])
ax6 = plt.subplot(gs[2, 2:4])
ax7 = plt.subplot(gs[3, 2:4])
    
#fig.suptitle(r'$I=$%s nA' % str(iamp),fontsize=20)

ax1.set_title(r'A',loc='left',fontsize=14)
ax2.set_title(r'B',loc='left',fontsize=14)
ax3.set_title(r'C',loc='left',fontsize=14)
ax4.set_title(r'D',loc='left',fontsize=14)
ax5.set_title(r'E',loc='left',fontsize=14)
ax6.set_title(r'F',loc='left',fontsize=14)
ax7.set_title(r'G',loc='left',fontsize=14)
ax1.set_title(r'Ca',fontsize=14)#loc='right',
ax2.set_title(r'SK',fontsize=14)#loc='right',
ax1.set_title(r'Trace, CaHVA only',fontsize=14)
ax2.set_title(r'Trace, CaHVA and SK',fontsize=14)
ax3.set_title(r'Voltage difference, CaHVA only',fontsize=14)
ax4.set_title(r'Voltage difference, CaHVA and SK',loc='right',fontsize=14)
ax5.set_title(r'$I_\mathregular{Ca}$, CaHVA only',fontsize=14)
ax6.set_title(r'$I_\mathregular{Ca}$, CaHVA and SK',fontsize=14)
ax7.set_title(r'$I_\mathregular{SK}$, CaHVA and SK',fontsize=14)
mycolors = ['tab:orange','tab:green','tab:red','tab:purple','tab:brown']
    
ax1.set_xlabel(r'$t$ (ms)',fontsize=12)
ax1.set_ylabel(r'$V$ (mV)',fontsize=12)
ax2.set_xlabel(r'$t$ (ms)',fontsize=12)
ax2.set_ylabel(r'$V$ (mV)',fontsize=12)
ax3.set_xlabel(r'$t$ (ms)',fontsize=12)
ax3.set_ylabel(r'$\Delta V$ (mV)',fontsize=12)
ax4.set_xlabel(r'$t$ (ms)',fontsize=12)
ax4.set_ylabel(r'$\Delta V$ (mV)',fontsize=12)
ax5.set_xlabel(r'$t$ (ms)',fontsize=12)
ax5.set_ylabel(r'$I_\mathregular{Ca}$ (nA)',fontsize=12)
ax6.set_xlabel(r'$t$ (ms)',fontsize=12)
ax6.set_ylabel(r'$I_\mathregular{Ca}$ (nA)',fontsize=12)
ax7.set_xlabel(r'$t$ (ms)',fontsize=12)
ax7.set_ylabel(r'$I_\mathregular{SK}$ (nA)',fontsize=12)

#t1_shift,V1_shifted_CaHVA

####################### CaHVA: ##########################################
ax1.axvline(x=peaktime_CaHVA,color='k',linestyle='--',linewidth=0.75)
ax1.plot(t5_shift, V5_shifted_CaHVA,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax1.plot(t4_shift, V4_shifted_CaHVA,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax1.plot(t3_shift, V3_shifted_CaHVA,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax1.plot(t2_shift, V2_shifted_CaHVA,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax1.plot(t1_shift, V1_shifted_CaHVA,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax1.set_xlim(left=t1_shift[0],right=t1_shift[-1])
ax1.legend(loc='upper right',ncol=1,fontsize=10)

ax3.axvline(x=peaktime_CaHVA,color='k',linestyle='--',linewidth=0.75)
ax3.plot(t1_shift, deltaV4,color=colors[0],label=r'$V_{%i\ \mathregular{mV}}-V_{%i\ \mathregular{mV}}$' % (Vshifts[0],Vshifts[2]))
ax3.plot(t1_shift, deltaV3,color=colors[1],label=r'$V_{%i\ \mathregular{mV}}-V_{%i\ \mathregular{mV}}$' % (Vshifts[1],Vshifts[2]))
ax3.plot(t1_shift, deltaV2,color=colors[3],label=r'$V_{%i\ \mathregular{mV}}-V_{%i\ \mathregular{mV}}$' % (Vshifts[3],Vshifts[2]))
ax3.plot(t1_shift, deltaV1,color=colors[4],label=r'$V_{%i\ \mathregular{mV}}-V_{%i\ \mathregular{mV}}$' % (Vshifts[4],Vshifts[2]))
ax3.set_xlim(left=peaktime_CaHVA,right=t1_shift[-1])
ax3.set_ylim(top=6,bottom=-7)
ax3.legend(loc='lower left',ncol=2,fontsize=10)

ax5.axvline(x=peaktime_CaHVA,color='k',linestyle='--',linewidth=0.75)
ax5.plot(t5_shift, I_Ca_HVA5_shifted_CaHVA,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax5.plot(t4_shift, I_Ca_HVA4_shifted_CaHVA,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax5.plot(t3_shift, I_Ca_HVA3_shifted_CaHVA,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax5.plot(t2_shift, I_Ca_HVA2_shifted_CaHVA,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax5.plot(t1_shift, I_Ca_HVA1_shifted_CaHVA,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax5.set_xlim(left=t1_shift[0],right=t1_shift[-1])
ax5.legend(loc='lower left',ncol=1,fontsize=10)


###
ax2.axvline(x=peaktime_CaHVA_SK,color='k',linestyle='--',linewidth=0.75)
ax2.plot(t5_shift_sk, V5_shifted_CaHVA_SK,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax2.plot(t4_shift_sk, V4_shifted_CaHVA_SK,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax2.plot(t3_shift_sk, V3_shifted_CaHVA_SK,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax2.plot(t2_shift_sk, V2_shifted_CaHVA_SK,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax2.plot(t1_shift_sk, V1_shifted_CaHVA_SK,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax2.set_xlim(left=t1_shift_sk[0],right=t1_shift_sk[-1])
ax2.legend(loc='upper right',ncol=1,fontsize=10)

ax4.axvline(x=peaktime_CaHVA_SK,color='k',linestyle='--',linewidth=0.75)
ax4.plot(t4_shift_sk, deltaV4_SK,color=colors[0],label=r'$V_{%i\ \mathregular{mV}}-V_{%i\ \mathregular{mV}}$' % (Vshifts[0],Vshifts[2]))
ax4.plot(t3_shift_sk, deltaV3_SK,color=colors[1],label=r'$V_{%i\ \mathregular{mV}}-V_{%i\ \mathregular{mV}}$' % (Vshifts[1],Vshifts[2]))
ax4.plot(t2_shift_sk, deltaV2_SK,color=colors[3],label=r'$V_{%i\ \mathregular{mV}}-V_{%i\ \mathregular{mV}}$' % (Vshifts[3],Vshifts[2]))
ax4.plot(t1_shift_sk, deltaV1_SK,color=colors[4],label=r'$V_{%i\ \mathregular{mV}}-V_{%i\ \mathregular{mV}}$' % (Vshifts[4],Vshifts[2]))
ax4.set_xlim(left=peaktime_CaHVA_SK,right=t1_shift_sk[-1])
ax4.set_ylim(top=6,bottom=-7)
ax4.legend(loc='lower left',ncol=2,fontsize=9)

ax6.axvline(x=peaktime_CaHVA_SK,color='k',linestyle='--',linewidth=0.75)
ax6.plot(t5_shift_sk, I_Ca_HVA_sk5_shifted_CaHVA_SK,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax6.plot(t4_shift_sk, I_Ca_HVA_sk4_shifted_CaHVA_SK,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax6.plot(t3_shift_sk, I_Ca_HVA_sk3_shifted_CaHVA_SK,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax6.plot(t2_shift_sk, I_Ca_HVA_sk2_shifted_CaHVA_SK,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax6.plot(t1_shift_sk, I_Ca_HVA_sk1_shifted_CaHVA_SK,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax6.set_xlim(left=t1_shift_sk[0],right=t1_shift_sk[-1])
ax6.legend(loc='lower left',ncol=1,fontsize=10)

#ax7.plot(t_nosk[istart:iend], I_SK_sk1_shifted_CaHVA_SK,label=r'0$\bar{g}_\mathregular{SK}$')
ax7.axvline(x=peaktime_CaHVA_SK,color='k',linestyle='--',linewidth=0.75)
ax7.plot(t5_shift_sk, I_SK_sk5_shifted_CaHVA_SK,color=colors[0],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[0])
ax7.plot(t4_shift_sk, I_SK_sk4_shifted_CaHVA_SK,color=colors[1],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[1])
ax7.plot(t3_shift_sk, I_SK_sk3_shifted_CaHVA_SK,color=colors[2],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[2])
ax7.plot(t2_shift_sk, I_SK_sk2_shifted_CaHVA_SK,color=colors[3],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[3])
ax7.plot(t1_shift_sk, I_SK_sk1_shifted_CaHVA_SK,color=colors[4],label=r'$V_\mathregular{shift}=$%i mV' % Vshifts[4])
ax7.set_xlim(left=t1_shift_sk[0],right=t1_shift_sk[-1])
ax7.legend(loc='upper left',ncol=1,fontsize=9)

fig.tight_layout()#(rect=[0, 0.03, 1, 0.95])
plt.savefig('Results/Soma10/Compare/trace_cainfo_idur'+str(idur)+'_iamp'+str(iamp)+'_CaHVA_SK_SHIFTS_aligned_smallset.png')
#plt.show()


################## Plotting total tracde and [Ca]in: ##################
fig = plt.figure(figsize=(10,8),dpi=300)#(figsize=(8,3),dpi=300)
    
gs = gridspec.GridSpec(2, 4)
    
ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[1, 0:2])
ax4 = plt.subplot(gs[1, 2:4])

ax1.set_title('A',loc='left')
ax2.set_title('B',loc='left')
ax3.set_title('C',loc='left')
ax4.set_title('D',loc='left')   

ax1.plot(t5, V5, label=r'$V_\mathregular{shift}$=%i mV' % Vshifts[0])
ax1.plot(t3, V3, label=r'$V_\mathregular{shift}$=%i mV' % Vshifts[2])
ax1.plot(t1, V1, label=r'$V_\mathregular{shift}$=%i mV' % Vshifts[4])
ax1.set_xlabel(r'$t$ (ms)')
ax1.set_ylabel(r'$V$ (mV)')
ax1.set_title('Trace, CaHVA only')
ax1.set_ylim(bottom=-85)
ax1.legend(loc='lower center',ncol=4,fontsize=7)

ax2.plot(t5, cai5, label=r'$V_\mathregular{shift}$=%i mV' % Vshifts[0])
ax2.plot(t3, cai3, label=r'$V_\mathregular{shift}$=%i mV' % Vshifts[2])
ax2.plot(t1, cai1, label=r'$V_\mathregular{shift}$=%i mV' % Vshifts[4])
ax2.set_xlabel(r'$t$ (ms)')
ax2.set_ylabel(r'$[\mathregular{Ca}^{2+}]_\mathregular{in}$ (mM)')
ax2.set_title(r'$[\mathregular{Ca}^{2+}]_\mathregular{in}$, CaHVA only')
ax2.legend(loc='lower right',ncol=2,fontsize=9)

ax3.plot(t_sk5, V_sk5, label=r'$V_\mathregular{shift}$=%i mV' % Vshifts[0])
ax3.plot(t_sk3, V_sk3, label=r'$V_\mathregular{shift}$=%i mV' % Vshifts[2])
ax3.plot(t_sk1, V_sk1, label=r'$V_\mathregular{shift}$=%i mV' % Vshifts[4])
ax3.set_xlabel(r'$t$ (ms)')
ax3.set_ylabel(r'$V$ (mV)')
ax3.set_title('Trace, CaHVA and SK')
ax3.set_ylim(bottom=-90)
ax3.legend(loc='lower center',ncol=4,fontsize=7)

ax4.plot(t_sk5, cai_sk5, label=r'$V_\mathregular{shift}$=%i mV' % Vshifts[0])
ax4.plot(t_sk3, cai_sk3, label=r'$V_\mathregular{shift}$=%i mV' % Vshifts[2])
ax4.plot(t_sk1, cai_sk1, label=r'$V_\mathregular{shift}$=%i mV' % Vshifts[4])
ax4.set_xlabel(r'$t$ (ms)')
ax4.set_ylabel(r'$[\mathregular{Ca}^{2+}]_\mathregular{in}$ (mM)')
ax4.set_title(r'$[\mathregular{Ca}^{2+}]_\mathregular{in}$, CaHVA and SK')
ax4.legend(loc='lower right',ncol=2,fontsize=9)

plt.tight_layout()
plt.savefig('Results/Soma10/Compare/trace_concentration_idur'+str(idur)+'_iamp'+str(iamp)+'_Ca_SK_SHIFTS_both.png')
#plt.show()

