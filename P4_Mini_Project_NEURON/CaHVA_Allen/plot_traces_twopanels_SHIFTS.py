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
    
fig = plt.figure(figsize=(10,4),dpi=300)#(figsize=(8,3),dpi=300)
    
gs = gridspec.GridSpec(1, 4)
    
ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
    
#fig.suptitle(r'Properties',fontsize=20)

ax1.set_title(r'A',loc='left',fontsize=14)
ax2.set_title(r'B',loc='left',fontsize=14)
ax1.set_title(r'Altering $\bar{g}_\mathregular{CaHVA}$, no SK',fontsize=14)#loc='right',
ax2.set_title(r'Altering $\bar{g}_\mathregular{SK}$, with CaHVA',fontsize=14)#loc='right',
    
ax1.set_xlabel(r'$t$ (ms)',fontsize=12)
ax1.set_ylabel(r'$V$ (mV)',fontsize=12)
ax2.set_xlabel(r'$t$ (ms)',fontsize=12)
ax2.set_ylabel(r'$V$ (mV)',fontsize=12)

idur = 1000
iamp = 0.003
gcahva = 0.2
gcahva_label = 1.0
Vshifts = [-20,-10,0,10,20]
NV = len(Vshifts)

ts = []
Vs = []

for i in range(NV):
    folder = 'Shift_gCaHVA_V/Vshift_'+str(Vshifts[i])+'/Results/Soma10/current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    filename = folder+'_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V.txt'
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

spdelta = 2200
epdelta = 2200
# Don't know these, should test out:
istart = int(5*N/8)+800
iend   = int(3*N/4)-7100
istart1 = istart+2800+spdelta
iend1   = iend-3541#+epdelta
istart2 = istart+1611+spdelta#+3000
iend2   = iend-4504+epdelta#-3400
istart3 = istart+536+spdelta#+2000
iend3   = iend-5579+epdelta#-4090
istart4 = istart+5675+spdelta#+8600
iend4   = iend-440+epdelta#+2500
istart5 = istart+1680+spdelta#+4200
iend5   = iend-4435+epdelta#-2829
ishift  = 2426#int(N/8)
ishift2 = 1237
ishift3 = 162 #-391
ishift4 = 5301 # Looking good
ishift5 = 1306

#'''
ax1.plot(t1[istart1:iend1], V1[istart1:iend1],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[0]))
ax1.plot(t2[istart2-ishift2:iend2-ishift2], V2[istart2:iend2],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[1]))
ax1.plot(t3[istart3-ishift3:iend3-ishift3], V3[istart3:iend3],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[2]))
ax1.plot(t4[istart4-ishift4:iend4-ishift4], V4[istart4:iend4],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[3]))
ax1.plot(t5[istart5-ishift5:iend5-ishift5], V5[istart5:iend5],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[4]))
ax1.axis([720.109,749.359,-78.2,51.1])
ax1.legend(loc='upper left',ncol=1,fontsize=11)

gsk      = 1.0

ts = []
Vs = []

for i in range(NV):
    folder = 'SK_Allen/Shift_gCaHVA_V/Vshift_'+str(Vshifts[i])+'/Results/Soma10/current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
    filename = folder+'_gSK'+str(gsk)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V.txt'
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

istart  = int(9*N/16)+1350
iend    = int(10*N/16)-2850
istart2 = int(istart*1.2)
iend2   = istart2+(iend-istart)
ishift1 = 3624
ishift2 = 1871
ishift3 = 2177
ishift4 = 5117
ishift5 = 2966 # -1915 is pretty good


print('t_sk1[istart]:',t_sk1[istart])
print('t_sk1[iend]:',t_sk1[iend])
print('t_sk1[istart2]:',t_sk1[istart2])
print('t_sk1[iend2]:',t_sk1[iend2])
print('V_sk3[istart2]:',V_sk3[istart2])
print('V_sk3[iend2]:',V_sk3[iend2])


ax2.plot(t_sk1[istart:iend], V_sk1[istart+ishift1:iend+ishift1],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[0]))
ax2.plot(t_sk2[istart:iend], V_sk2[istart+ishift2:iend+ishift2],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[1]))
ax2.plot(t_sk3[istart:iend], V_sk3[istart+ishift3:iend+ishift3],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[2]))
ax2.plot(t_sk4[istart:iend], V_sk4[istart+ishift4:iend+ishift4],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[3]))
ax2.plot(t_sk5[istart:iend], V_sk5[istart2+ishift5:iend2+ishift5],label=r'$V_\mathregular{shift}=%s$ mV' % str(Vshifts[4]))
ax2.axis([634.922,671.476,-78.2,51.1])
ax2.legend(loc='upper left',ncol=1,fontsize=11)

plt.tight_layout()
plt.savefig('Results/Soma10/Compare/trace_twopanels_idur'+str(idur)+'_iamp'+str(iamp)+'_aligned_SHIFTS.png')
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
    plt.savefig('Results/Soma10/Compare/trace_twopanels_idur'+str(idur)+'_iamp'+str(iamp)+'_unaligned_SHIFTS.png')
    plt.show()

