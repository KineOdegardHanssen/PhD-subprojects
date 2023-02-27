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
gsk      = 1.0
gsk2     = 1.1
gcahva   = 0.2
gcahva2  = 1.0
gcahva_plot  = 1.0
gcahva_plot2 = 5.0

folderbase    = 'Results/Soma10/current_idur'+str(idur)+'_iamp'+str(iamp)+'/'
filename_nosk = folderbase+'_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V.txt'
filename_sk_1 = 'SK_Allen/'+folderbase+'_gSK'+str(gsk)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V.txt'
filename_sk_2 = 'SK_Allen/'+folderbase+'_gSK'+str(gsk2)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V.txt'

t_sk1, V_sk1 = unpack(filename_sk_1)
t_sk2, V_sk2 = unpack(filename_sk_2)
t_nosk, V_nosk   = unpack(filename_nosk)

N = len(t_sk1)

istart = int(5*N/8)+1300
iend   = int(3*N/4)-2800
ishift = 2426#int(N/8)
ishift2 = 766 #OK: 765

plt.figure(figsize=(6,5))
plt.plot(t_nosk[istart:iend], V_nosk[istart+ishift2:iend+ishift2],label=r'0$\bar{g}_\mathregular{SK}$, %s$\bar{g}_\mathregular{CaHVA}$' % str(gcahva_plot))
plt.plot(t_sk1[istart:iend], V_sk1[istart+ishift:iend+ishift],label=r'%s$\bar{g}_\mathregular{SK}$, %s$\bar{g}_\mathregular{CaHVA}$' % (str(gsk),str(gcahva_plot)))
plt.plot(t_sk2[istart:iend], V_sk2[istart+ishift:iend+ishift],label=r'%s$\bar{g}_\mathregular{SK}$, %s$\bar{g}_\mathregular{CaHVA}$' % (str(gsk2),str(gcahva_plot)))
plt.xlabel('t (ms)')
plt.ylabel('V (mV)')
plt.legend(loc='upper left',ncol=1,fontsize=11)
plt.tight_layout()
plt.savefig('Results/Soma10/Compare/trace_withandwithoutSK_idur'+str(idur)+'_iamp'+str(iamp)+'_gcahva_%.1f_gsk_%.1f_updated.png' % (gcahva_plot,gsk))
plt.show()

plt.figure(figsize=(6,5))
plt.plot(t_nosk, V_nosk,label=r'0$\bar{g}_\mathregular{SK}$, %s$\bar{g}_\mathregular{CaHVA}$' % str(gcahva_plot))
plt.plot(t_sk1, V_sk1,label=r'%s$\bar{g}_\mathregular{SK}$, %s$\bar{g}_\mathregular{CaHVA}$' % (str(gsk),str(gcahva_plot)))
plt.plot(t_sk2, V_sk2,label=r'%s$\bar{g}_\mathregular{SK}$, %s$\bar{g}_\mathregular{CaHVA}$' % (str(gsk2),str(gcahva_plot)))
plt.xlabel('t (ms)')
plt.ylabel('V (mV)')
plt.legend(loc='upper left',ncol=1,fontsize=11)
plt.tight_layout()
plt.savefig('Results/Soma10/Compare/trace_withandwithoutSK_idur'+str(idur)+'_iamp'+str(iamp)+'_gcahva_%.1f_gsk_%.1f_unaligned_updated.png' % (gcahva_plot,gsk))
plt.show()

