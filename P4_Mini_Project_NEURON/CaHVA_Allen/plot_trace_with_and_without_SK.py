import matplotlib.pyplot as plt
import numpy as np

def unpack(filename):
    data = np.loadtxt(filename)

    # Time is the first column
    x = data[:, 0]
    # Voltage is the second column
    y = data[:, 1]
    return x, y


idur = 1000
iamp = 0.002
gsk     = 1.0
gcahva1 = 1.0
gcahva2 = 1.0

folderbase    = 'Results/Soma10/current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
filename_nosk = folderbase+'_gCaHVA'+str(gcahva1)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V.txt'
filename_sk   = 'SK_Allen/'+folderbase+'_gSK'+str(gsk)+'p_gCaHVA'+str(gcahva2)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V.txt'

t_sk, V_sk = unpack(filename_sk)
t_nosk, V_nosk = unpack(filename_nosk)

plt.figure(figsize=(6,5))
plt.plot(t_nosk, V_nosk,label=r'0$\bar{g}_\mathregular{SK}$')
plt.plot(t_sk, V_sk,label=r'1$\bar{g}_\mathregular{SK}$')
plt.xlabel('t (ms)')
plt.ylabel('V (mV)')
plt.legend(loc='upper left',ncol=1,fontsize=11)
plt.tight_layout()
plt.savefig('Results/Soma10/Compare/trace_withandwithoutSK_idur'+str(idur)+'_iamp'+str(iamp)+'_gcahva1_%.1f_gcahva2_%.1f_gsk_%.1f.png' % (gcahva1,gcahva2,gsk))


