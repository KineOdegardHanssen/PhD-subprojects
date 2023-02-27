import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np

def unpack(filename):
    data = np.loadtxt(filename)

    # Time is the first column
    x = data[:, 0]
    # Voltage is the second column
    y = data[:, 1]
    return x, y

plotitall = False # True # 


plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
    
fig = plt.figure(figsize=(10,8),dpi=300)#(figsize=(8,3),dpi=300)
    
gs = gridspec.GridSpec(2, 4)

ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[1, 0:2])
ax4 = plt.subplot(gs[1, 2:4])
    
#fig.suptitle(r'Properties',fontsize=20)

ax1.set_title(r'A',loc='left',fontsize=14)
ax2.set_title(r'B',loc='left',fontsize=14)
ax3.set_title(r'C',loc='left',fontsize=14)
ax4.set_title(r'D',loc='left',fontsize=14)
ax1.set_title(r'CaHVA, no SK',fontsize=14)#loc='right',
ax2.set_title(r'CaHVA and SK',fontsize=14)#loc='right',
ax3.set_title(r'CaHVA, no SK',fontsize=14)#loc='right',
ax4.set_title(r'CaHVA and SK',fontsize=14)#loc='right',
    
ax1.set_xlabel(r'$t$ (ms)',fontsize=12)
ax1.set_ylabel(r'$V$ (mV)',fontsize=12)
ax2.set_xlabel(r'$t$ (ms)',fontsize=12)
ax2.set_ylabel(r'$V$ (mV)',fontsize=12)
ax3.set_xlabel(r'$t$ (ms)',fontsize=12)
ax3.set_ylabel(r'$V_{\mathregular{shift}=i}$-$V_{\mathregular{shift}=-20\ \mathregular{mV}}$',fontsize=16)
ax4.set_xlabel(r'$t$ (ms)',fontsize=12)
ax4.set_ylabel(r'$V_{\mathregular{shift}=i}$-$V_{\mathregular{shift}=-20\ \mathregular{mV}}$',fontsize=16)
mycolors = ['tab:orange','tab:green','tab:red','tab:purple']

idur = 1000
iamp = 0.02
gcahva = 0.2
gcahva_label = 1.0
Vshifts = [-20,-10,0,10,20]
NV = len(Vshifts)

ts = []
Vs = []

for i in range(NV):
    folder = 'Shift_gCaHVA_V/Vshift_'+str(Vshifts[i])+'/Results/Soma10/current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    filename = folder+'_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V_cainfo.txt'
    t, V = unpack(filename)
    N = len(t)
    #ax1.plot(t[int(N/2):], V[int(N/2):],label=r'%s$\bar{g}$' % gcahvas_label[i])
    ts.append(t)
    Vs.append(V)
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

shiftall = -100#200
spdelta = 2200
epdelta = 2200
# Don't know these, should test out:
istart = int(9*N/16)+938+shiftall
iend   = int(5*N/8)-848+shiftall
istart1 = istart+2642
iend1   = iend-3473
istart2 = istart+1511+spdelta#+3000
iend2   = iend-4604+epdelta#-3400
istart3 = istart+577+spdelta#+2000
iend3   = iend-5538+epdelta#-4090
istart4 = istart+6019+spdelta#+8600
iend4   = iend-96+epdelta#+2500
istart5 = istart+1858+spdelta#+4200
iend5   = iend-4257+epdelta#-2829
ishift1 = 0#-200#2426
ishift2 = 1069
ishift3 = 135
ishift4 = 5577#5575
ishift5 = 1416

t1_shift = t1[istart1:iend1]
V1_CaHVA_shifted = V1[istart1:iend1]
V2_CaHVA_shifted = V2[istart2:iend2]
V3_CaHVA_shifted = V3[istart3:iend3]
V4_CaHVA_shifted = V4[istart4:iend4]
V5_CaHVA_shifted = V5[istart5:iend5]

ax1.plot(t1[istart1:iend1], V1[istart1:iend1],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[0]))
ax1.plot(t2[istart2-ishift2:iend2-ishift2], V2[istart2:iend2],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[1]))
ax1.plot(t3[istart3-ishift3:iend3-ishift3], V3[istart3:iend3],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[2]))
ax1.plot(t4[istart4-ishift4:iend4-ishift4], V4[istart4:iend4],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[3]))
ax1.plot(t5[istart5-ishift5:iend5-ishift5], V5[istart5:iend5],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[4]))
ax1.set_xlim(left=t1[istart1],right=t1[iend1])
#ax1.axis([720.109,749.359,-78.2,51.1])
ax1.legend(loc='upper right',ncol=1,fontsize=11)

#'''
deltaV1 = V2_CaHVA_shifted-V1_CaHVA_shifted
deltaV2 = V3_CaHVA_shifted-V1_CaHVA_shifted
deltaV3 = V4_CaHVA_shifted-V1_CaHVA_shifted
deltaV4 = V5_CaHVA_shifted-V1_CaHVA_shifted


print('len(t1_shift):',len(t1_shift))
print('len(deltaV1):',len(deltaV1))
print('len(deltaV2):',len(deltaV2))
print('len(deltaV3):',len(deltaV3))
print('len(deltaV4):',len(deltaV4))

ax3.plot(t1_shift, deltaV1,label=r'$i=%s$ mV' % str(Vshifts[1]),color=mycolors[0])
ax3.plot(t1_shift, deltaV2,label=r'$i=%s$ mV' % str(Vshifts[2]),color=mycolors[1])
ax3.plot(t1_shift, deltaV3,label=r'$i=%s$ mV' % str(Vshifts[3]),color=mycolors[2])
ax3.plot(t1_shift, deltaV4,label=r'$i=%s$ mV' % str(Vshifts[4]),color=mycolors[3])
ax3.set_xlim(left=t1[istart1],right=t1[iend1])
#ax3.axis([720.109,749.359,-78.2,51.1])
ax3.legend(loc='lower left',ncol=1,fontsize=11)
#'''

gsk      = 1.0

ts = []
Vs = []

for i in range(NV):
    folder = 'SK_Allen/Shift_gCaHVA_V/Vshift_'+str(Vshifts[i])+'/Results/Soma10/current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    filename = folder+'_gSK'+str(gsk)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V_writecainfo.txt'
    t, V = unpack(filename)
    N = len(t)
    #ax1.plot(t[int(N/2):], V[int(N/2):],label=r'%s$\bar{g}$' % gcahvas_label[i])
    ts.append(t)
    Vs.append(V)

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

Nt = len(t_sk1)

shorten1 = 0#467
istart = int(9*N/16)+800
iend   = int(5*N/8)-7100
istart2 = int(istart*1.2)
iend2   = istart2+(iend-istart)
istart  = istart+shorten1
iend    = iend-shorten1-1
istart2 = istart2+shorten1
iend2   = iend2-shorten1-1
ishift1 = 3784
ishift2 = 1683
ishift3 = 2018
ishift4 = 4401
ishift5 = 3109 # -1915 is pretty good


print('t_sk1[istart]:',t_sk1[istart])
print('t_sk1[iend]:',t_sk1[iend])
print('t_sk1[istart2]:',t_sk1[istart2])
print('t_sk1[iend2]:',t_sk1[iend2])
print('V_sk3[istart2]:',V_sk3[istart2])
print('V_sk3[iend2]:',V_sk3[iend2])

t1_sk_shift = t_sk1[istart:iend]
V1_CaHVA_SK_shifted = V_sk1[istart+ishift1:iend+ishift1]
V2_CaHVA_SK_shifted = V_sk2[istart+ishift2:iend+ishift2]
V3_CaHVA_SK_shifted = V_sk3[istart+ishift3:iend+ishift3]
V4_CaHVA_SK_shifted = V_sk4[istart+ishift4:iend+ishift4]
V5_CaHVA_SK_shifted = V_sk5[istart2+ishift5:iend2+ishift5]

deltaV1_sk = V2_CaHVA_SK_shifted-V1_CaHVA_SK_shifted
deltaV2_sk = V3_CaHVA_SK_shifted-V1_CaHVA_SK_shifted
deltaV3_sk = V4_CaHVA_SK_shifted-V1_CaHVA_SK_shifted
deltaV4_sk = V5_CaHVA_SK_shifted-V1_CaHVA_SK_shifted

'''
print('len(V1_CaHVA_SK_shifted):',len(V1_CaHVA_SK_shifted))
print('len(V2_CaHVA_SK_shifted):',len(V2_CaHVA_SK_shifted))
print('len(V3_CaHVA_SK_shifted):',len(V3_CaHVA_SK_shifted))
print('len(V4_CaHVA_SK_shifted):',len(V4_CaHVA_SK_shifted))
print('len(V5_CaHVA_SK_shifted):',len(V5_CaHVA_SK_shifted))
'''

print('len(t1_sk_shift):',len(t1_sk_shift))
print('len(deltaV1_sk):',len(deltaV1_sk))
print('len(deltaV2_sk):',len(deltaV2_sk))
print('len(deltaV3_sk):',len(deltaV3_sk))
print('len(deltaV4_sk):',len(deltaV4_sk))

'''
print('deltaV1_sk:',deltaV1_sk)
print('deltaV2_sk:',deltaV2_sk)
print('deltaV3_sk:',deltaV3_sk)
print('deltaV4_sk:',deltaV4_sk)
'''

ax2.plot(t_sk1[istart:iend], V_sk1[istart+ishift1:iend+ishift1],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[0]))
ax2.plot(t_sk2[istart:iend], V_sk2[istart+ishift2:iend+ishift2],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[1]))
ax2.plot(t_sk3[istart:iend], V_sk3[istart+ishift3:iend+ishift3],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[2]))
ax2.plot(t_sk4[istart:iend], V_sk4[istart+ishift4:iend+ishift4],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[3]))
ax2.plot(t_sk5[istart:iend], V_sk5[istart2+ishift5:iend2+ishift5],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[4]))
ax2.set_xlim(left=t_sk1[istart],right=t_sk1[iend])
#ax2.axis([634.922,671.476,-78.2,51.1])
ax2.legend(loc='upper left',ncol=1,fontsize=11)



ax4.plot(t1_sk_shift, deltaV1_sk,label=r'$i=%s$ mV' % str(Vshifts[1]),color=mycolors[0])
ax4.plot(t1_sk_shift, deltaV2_sk,label=r'$i=%s$ mV' % str(Vshifts[2]),color=mycolors[1])
ax4.plot(t1_sk_shift, deltaV3_sk,label=r'$i=%s$ mV' % str(Vshifts[3]),color=mycolors[2])
ax4.plot(t1_sk_shift, deltaV4_sk,label=r'$i=%s$ mV' % str(Vshifts[4]),color=mycolors[3])
ax4.set_xlim(left=t1_sk_shift[0],right=t1_sk_shift[-1])
#ax3.axis([720.109,749.359,-78.2,51.1])
ax4.legend(loc='lower left',ncol=1,fontsize=11)


plt.tight_layout()
plt.savefig('Results/Soma10/Compare/trace_fourpanels_idur'+str(idur)+'_iamp'+str(iamp)+'_aligned_SHIFTS.png')
plt.show()

if plotitall==True:
    fig = plt.figure(figsize=(12,6),dpi=300)
     
    gs = gridspec.GridSpec(1, 4)
    
    ax1 = plt.subplot(gs[0, 0:2])
    ax2 = plt.subplot(gs[0, 2:4])
    
    ax1.set_title(r'A',loc='left',fontsize=14)
    ax2.set_title(r'B',loc='left',fontsize=14)
    ax1.set_title(r'Altering $\bar{g}_\mathregular{CaHVA}$, no SK',fontsize=14)
    ax2.set_title(r'Altering $\bar{g}_\mathregular{SK}$',fontsize=14)
    
    ax1.set_xlabel(r'$t$ (ms)',fontsize=12)
    ax1.set_ylabel(r'$V$ (mV)',fontsize=12)
    ax2.set_xlabel(r'$t$ (ms)',fontsize=12)
    ax2.set_ylabel(r'$V$ (mV)',fontsize=12)
    
    ax1.plot(t1, V1,label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[0]))
    ax1.plot(t2, V2,label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[1]))
    ax1.plot(t3, V3,label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[2]))
    ax1.plot(t4, V4,label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[3]))
    ax1.plot(t5, V5,label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[4]))
    ax1.legend(loc='upper left',ncol=1,fontsize=11)
    
    ax2.plot(t_sk1, V_sk1,label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[0]))
    ax2.plot(t_sk2, V_sk2,label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[1]))
    ax2.plot(t_sk3, V_sk3,label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[2]))
    ax2.plot(t_sk4, V_sk4,label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[3]))
    ax2.plot(t_sk5, V_sk5,label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[4]))
    ax2.legend(loc='upper left',ncol=1,fontsize=11)
    plt.tight_layout()
    plt.savefig('Results/Soma10/Compare/trace_fourpanels_idur'+str(idur)+'_iamp'+str(iamp)+'_unaligned_SHIFTS.png')
    plt.show()

