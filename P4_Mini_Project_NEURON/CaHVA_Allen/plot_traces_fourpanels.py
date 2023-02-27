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
ax3.set_title(r'A',loc='left',fontsize=14)
ax4.set_title(r'B',loc='left',fontsize=14)
ax1.set_title(r'Altering $\bar{g}_\mathregular{CaHVA}$, no SK',fontsize=14)#loc='right',
ax2.set_title(r'Altering $\bar{g}_\mathregular{SK}$, with CaHVA',fontsize=14)#loc='right',
ax3.set_title(r'Altering $\bar{g}_\mathregular{CaHVA}$, no SK',fontsize=14)#loc='right',
ax4.set_title(r'Altering $\bar{g}_\mathregular{SK}$, with CaHVA',fontsize=14)#loc='right',
mycolors = ['tab:orange','tab:green','tab:red','tab:purple','tab:brown']
    
ax1.set_xlabel(r'$t$ (ms)',fontsize=12)
ax1.set_ylabel(r'$V$ (mV)',fontsize=12)
ax2.set_xlabel(r'$t$ (ms)',fontsize=12)
ax2.set_ylabel(r'$V$ (mV)',fontsize=12)
ax3.set_xlabel(r'$t$ (ms)',fontsize=12)
ax3.set_ylabel(r'$\Delta V$ (mV)',fontsize=12)
ax4.set_xlabel(r'$t$ (ms)',fontsize=12)
ax4.set_ylabel(r'$\Delta V$ (mV)',fontsize=12)

idur = 1000
iamp = 0.003
gcahvas = [0.2,0.5,1.0,1.5,2.0]
gcahvas_label = [1.0,2.5,5.0,7.5,10.0]
Ng = len(gcahvas)

folderbase    = 'Results/Soma10/current_idur'+str(idur)+'_iamp'+str(iamp)+'/'

ts = []
Vs = []

for i in range(Ng):
    filename = folderbase+'_gCaHVA'+str(gcahvas[i])+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V.txt'
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

istart = int(5*N/8)+800
iend   = int(3*N/4)-7100
istart1 = istart+3259#
iend1   = iend-3898
istart2 = istart+3757#
iend2   = iend-3400
istart3 = istart+2968#
iend3   = iend-4189
istart4 = istart+9200#
iend4   = iend+2043
istart5 = istart+3874#
iend5   = iend-3283
ishift  = 2426#int(N/8)
ishift2 = 498
ishift3 = -291
ishift4 = 5941 # Looking good
ishift5 = 615

t1_shift = t1[istart1:iend1]
V1_shifted_CaHVA = V1[istart1:iend1]
V2_shifted_CaHVA = V2[istart2:iend2]
V3_shifted_CaHVA = V3[istart3:iend3]
V4_shifted_CaHVA = V4[istart4:iend4]
V5_shifted_CaHVA = V5[istart5:iend5]

deltaV1 = V2_shifted_CaHVA-V1_shifted_CaHVA
deltaV2 = V3_shifted_CaHVA-V1_shifted_CaHVA
deltaV3 = V4_shifted_CaHVA-V1_shifted_CaHVA
deltaV4 = V5_shifted_CaHVA-V1_shifted_CaHVA

print('len(V1_shifted_CaHVA):',len(V1_shifted_CaHVA))
print('len(V2_shifted_CaHVA):',len(V2_shifted_CaHVA))
print('len(V3_shifted_CaHVA):',len(V3_shifted_CaHVA))
print('len(V4_shifted_CaHVA):',len(V4_shifted_CaHVA))
print('len(V5_shifted_CaHVA):',len(V5_shifted_CaHVA))

#'''
ax1.plot(t1[istart1:iend1], V1[istart1:iend1],label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[0]))
ax1.plot(t2[istart2-ishift2:iend2-ishift2], V2[istart2:iend2],label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[1]))
ax1.plot(t3[istart3-ishift3:iend3-ishift3], V3[istart3:iend3],label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[2]))
ax1.plot(t4[istart4-ishift4:iend4-ishift4], V4[istart4:iend4],label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[3]))
ax1.plot(t5[istart5-ishift5:iend5-ishift5], V5[istart5:iend5],label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[4]))
ax1.set_xlim(left=t1[istart1],right=t1[iend1])
#ax1.axis([723,746,-78.2,51.1])
ax1.legend(loc='upper left',ncol=1,fontsize=11)

ax3.plot(t1_shift, deltaV1,label=r'$V_{%s\bar{g}_\mathregular{CaHVA}}$-$V_{%s\bar{g}_\mathregular{CaHVA}}$' % (str(gcahvas_label[1]),str(gcahvas_label[0])),color=mycolors[0])
ax3.plot(t1_shift, deltaV2,label=r'$V_{%s\bar{g}_\mathregular{CaHVA}}$-$V_{%s\bar{g}_\mathregular{CaHVA}}$' % (str(gcahvas_label[2]),str(gcahvas_label[0])),color=mycolors[1])
ax3.plot(t1_shift, deltaV3,label=r'$V_{%s\bar{g}_\mathregular{CaHVA}}$-$V_{%s\bar{g}_\mathregular{CaHVA}}$' % (str(gcahvas_label[3]),str(gcahvas_label[0])),color=mycolors[2])
ax3.plot(t1_shift, deltaV4,label=r'$V_{%s\bar{g}_\mathregular{CaHVA}}$-$V_{%s\bar{g}_\mathregular{CaHVA}}$' % (str(gcahvas_label[4]),str(gcahvas_label[0])),color=mycolors[3])
ax3.set_xlim(left=t1[istart1],right=t1[iend1])
#ax3.axis([720.109,749.359,-78.2,51.1])
ax3.legend(loc='lower left',ncol=1,fontsize=11)


gsk      = 1.0
gsk2     = 2.0
gsk3     = 3.0
gsk4     = 5.0
gsk5     = 10.0
gcahva   = 0.2
gcahva2  = 1.0
gcahva_plot  = 1.0
gcahva_plot2 = 5.0

folderbase    = 'Results/Soma10/current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
filename_nosk = folderbase+'_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V.txt'
filename_sk_1 = 'SK_Allen/'+folderbase+'_gSK'+str(gsk)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V.txt'
filename_sk_2 = 'SK_Allen/'+folderbase+'_gSK'+str(gsk2)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V.txt'
filename_sk_3 = 'SK_Allen/'+folderbase+'_gSK'+str(gsk3)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V.txt'
filename_sk_4 = 'SK_Allen/'+folderbase+'_gSK'+str(gsk4)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V.txt'
filename_sk_5 = 'SK_Allen/'+folderbase+'_gSK'+str(gsk5)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V.txt'

t_sk1, V_sk1 = unpack(filename_sk_1)
t_sk2, V_sk2 = unpack(filename_sk_2)
t_sk3, V_sk3 = unpack(filename_sk_3)
t_sk4, V_sk4 = unpack(filename_sk_4)
t_sk5, V_sk5 = unpack(filename_sk_5)
t_nosk, V_nosk   = unpack(filename_nosk)

N = len(t_sk1)

shorten1 = 989
istart  = int(9*N/16)+1350
iend    = int(10*N/16)-2850
istart2 = int(istart*1.2)
iend2   = istart2+(iend-istart)
istart  = istart+shorten1-1
iend    = iend-shorten1
istart2 = istart2+shorten1-1
iend2   = iend2-shorten1
ishift  = 2413
ishift2 = 766
ishift3 = 1671
ishift4 = 4009
ishift5 = 3158
ishift6 = -1916 # -1915 is pretty good


print('t[istart]:',t_nosk[istart])
print('t[iend]:',t_nosk[iend])
print('t[istart2]:',t_nosk[istart2])
print('t[iend2]:',t_nosk[iend2])
print('V_sk5[istart2]:',V_sk3[istart2])
print('V_sk3[iend2]:',V_sk3[iend2])


t1_shift_sk = t_nosk[istart:iend]
V1_shifted_CaHVA_SK = V_nosk[istart+ishift2:iend+ishift2]
V2_shifted_CaHVA_SK = V_sk1[istart+ishift:iend+ishift]
V3_shifted_CaHVA_SK = V_sk2[istart+ishift3:iend+ishift3]
V4_shifted_CaHVA_SK = V_sk3[istart+ishift4:iend+ishift4]
V5_shifted_CaHVA_SK = V_sk4[istart+ishift5:iend+ishift5]
V6_shifted_CaHVA_SK = V_sk5[istart2+ishift6:iend2+ishift6]

print('len(V1_shifted_CaHVA_SK):',len(V1_shifted_CaHVA_SK))

deltaV1 = V2_shifted_CaHVA_SK-V1_shifted_CaHVA_SK
deltaV2 = V3_shifted_CaHVA_SK-V1_shifted_CaHVA_SK
deltaV3 = V4_shifted_CaHVA_SK-V1_shifted_CaHVA_SK
deltaV4 = V5_shifted_CaHVA_SK-V1_shifted_CaHVA_SK
deltaV5 = V6_shifted_CaHVA_SK-V1_shifted_CaHVA_SK


ax2.plot(t_nosk[istart:iend], V_nosk[istart+ishift2:iend+ishift2],label=r'0$\bar{g}_\mathregular{SK}$')
ax2.plot(t_sk1[istart:iend], V_sk1[istart+ishift:iend+ishift],label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk))
ax2.plot(t_sk2[istart:iend], V_sk2[istart+ishift3:iend+ishift3],label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk2))
ax2.plot(t_sk3[istart:iend], V_sk3[istart+ishift4:iend+ishift4],label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk3))
ax2.plot(t_sk4[istart:iend], V_sk4[istart+ishift5:iend+ishift5],label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk4))
ax2.plot(t_sk5[istart:iend], V_sk5[istart2+ishift6:iend2+ishift6],label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk5))
ax2.set_xlim(left=t_nosk[istart],right=t_nosk[iend])
#ax2.axis([634.84,671.6,-78.2,51.1])
ax2.legend(loc='upper left',ncol=1,fontsize=11)


ax4.plot(t1_shift_sk, deltaV1,label=r'$V_{%s\bar{g}_\mathregular{SK}}$-$V_{0\bar{g}_\mathregular{SK}}$' % str(gsk),color=mycolors[0])
ax4.plot(t1_shift_sk, deltaV2,label=r'$V_{%s\bar{g}_\mathregular{SK}}$-$V_{0\bar{g}_\mathregular{SK}}$' % str(gsk2),color=mycolors[1])
ax4.plot(t1_shift_sk, deltaV3,label=r'$V_{%s\bar{g}_\mathregular{SK}}$-$V_{0\bar{g}_\mathregular{SK}}$' % str(gsk3),color=mycolors[2])
ax4.plot(t1_shift_sk, deltaV4,label=r'$V_{%s\bar{g}_\mathregular{SK}}$-$V_{0\bar{g}_\mathregular{SK}}$' % str(gsk4),color=mycolors[3])
ax4.plot(t1_shift_sk, deltaV5,label=r'$V_{%s\bar{g}_\mathregular{SK}}$-$V_{0\bar{g}_\mathregular{SK}}$' % str(gsk5),color=mycolors[4])
ax4.set_xlim(left=t1_shift_sk[0],right=t1_shift_sk[-1])
#ax3.axis([720.109,749.359,-78.2,51.1])
ax4.legend(loc='upper right',ncol=1,fontsize=11)


print('len(t1_shift):',len(t1_shift))
print('len(t1_shift_sk):',len(t1_shift_sk))

plt.tight_layout()
plt.savefig('Results/Soma10/Compare/trace_fourpanels_idur'+str(idur)+'_iamp'+str(iamp)+'_aligned.png')
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
    
    ax1.plot(t1, V1,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[0]))
    ax1.plot(t2, V2,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[1]))
    ax1.plot(t3, V3,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[2]))
    ax1.plot(t4, V4,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[3]))
    ax1.plot(t5, V5,label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[4]))
    ax1.legend(loc='upper left',ncol=1,fontsize=11)

    ax2.plot(t_nosk, V_nosk,label=r'0$\bar{g}_\mathregular{SK}$' )
    ax2.plot(t_sk1, V_sk1,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk))
    ax2.plot(t_sk2, V_sk2,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk2))
    ax2.plot(t_sk3, V_sk3,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk3))
    ax2.plot(t_sk4, V_sk4,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk4))
    ax2.plot(t_sk5, V_sk5,label=r'%s$\bar{g}_\mathregular{SK}$' % str(gsk5))
    ax2.legend(loc='upper left',ncol=1,fontsize=11)
    plt.tight_layout()
    plt.savefig('Results/Soma10/Compare/trace_fourpanels_idur'+str(idur)+'_iamp'+str(iamp)+'_unaligned.png')
    plt.show()

