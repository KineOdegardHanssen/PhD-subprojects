import matplotlib.pyplot as plt
import numpy as np

def unpack(filename):
    data = np.loadtxt(filename)

    # Time is the first column
    x = data[:, 0]
    # Voltage is the second column
    y = data[:, 1]
    return x, y


idur = 300
iamp = 0.002
gcahvas = [0.2,0.5,1.0,1.5,2.0]
gcahvas_label = [1.0,2.5,5.0,7.5,10.0]
Ng = len(gcahvas)

folderbase    = 'Results/Soma10/current_idur'+str(idur)+'_iamp'+str(iamp)+'/'

ts = []
Vs = []

plt.figure(figsize=(6,5))
for i in range(Ng):
    filename = folderbase+'_gCaHVA'+str(gcahvas[i])+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V.txt'
    t, V = unpack(filename)
    N = len(t)
    plt.plot(t[int(N/2):], V[int(N/2):],label=r'%s$\bar{g}$' % gcahvas_label[i])
    ts.append(t)
    Vs.append(V)
plt.xlabel('t (ms)')
plt.ylabel('V (mV)')
plt.legend(loc='upper left',ncol=1,fontsize=11)
plt.tight_layout()
plt.savefig('Results/Soma10/Compare/trace_withandwithoutSK_idur'+str(idur)+'_iamp'+str(iamp)+'_varygcahva_unaligned.png')
plt.show()


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

istart = int(5*N/8)#+1300
iend   = int(3*N/4)#-2800
istart1 = int(5*N/8)+2800
iend1   = int(3*N/4)-3315
istart2 = int(5*N/8)+2600
iend2   = int(3*N/4)-3500
istart3 = int(5*N/8)+2000
iend3   = int(3*N/4)-4090
istart4 = int(5*N/8)+8600
iend4   = int(3*N/4)+2500
istart5 = int(5*N/8)+3300
iend5   = int(3*N/4)-2829
ishift = 2426#int(N/8)
ishift2 = -179
ishift3 = -767
ishift4 = 5828
ishift5 = 495

#'''
plt.figure(figsize=(6,5))
plt.plot(t1[istart1:iend1], V1[istart1:iend1],label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[0]))
plt.plot(t2[istart2-ishift2:iend2-ishift2], V2[istart2:iend2],label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[1]))
plt.plot(t3[istart3-ishift3:iend3-ishift3], V3[istart3:iend3],label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[2]))
plt.plot(t4[istart4-ishift4:iend4-ishift4], V4[istart4:iend4],label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[3]))
plt.plot(t5[istart5-ishift5:iend5-ishift5], V5[istart5:iend5],label=r'%s$\bar{g}_\mathregular{CaHVA}$' % str(gcahvas_label[4]))
plt.xlabel('t (ms)')
plt.ylabel('V (mV)')
plt.axis([278.625,281.472,-72.2,49.5])
plt.legend(loc='upper left',ncol=1,fontsize=11)
plt.tight_layout()
plt.savefig('Results/Soma10/Compare/trace_withandwithoutSK_idur'+str(idur)+'_iamp'+str(iamp)+'_varygcahvas_aligned.png')
plt.show()
#'''



