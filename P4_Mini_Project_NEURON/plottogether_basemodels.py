import numpy as np
import matplotlib.pyplot as plt

def unpack(filename):
    data = np.loadtxt(filename)

    # Time is the first column
    x = data[:, 0]
    # Voltage is the second column
    y = data[:, 1]
    return x, y

figname = 'compare_Allenbase_Hjorthbase_I0.01nA.png'

# Trace file
allenfile = 'Allen_NaV_Kv2like/Results/Soma10/current_idur100_iamp0.01/Epasshift-20_gNaV0.5p_gKv2like0.5p_somaonly_cm1.0_idur100_iamp0.01_dtexp-7_V.txt'
hjorthfile = 'Hjorth2020/FS_Na_K/Results/Soma10/current_idur100_iamp0.01/somaonly_cm1.0_idur100_iamp0.01_dtexp-7_vinit-86_trec-600.0_V.txt'

t_allen, V_allen = unpack(allenfile)
t_hjorth, V_hjorth = unpack(hjorthfile)

# fI-curve
fI_hjorth_file = 'Hjorth2020/FS_Na_K/Results/Soma10/somaHjorth_idur1000_varyiamp_manual_cm1.0_Nspikes_vs_I_s500.txt'
fI_allen_file  = 'Allen_NaV_Kv2like/Results/Soma10/somaAllen_idur1000_varyiamp_manual_cm1.0_Epasshift-20_gNaV0.5p_gKv2like0.5p_Nspikes_vs_I_s500.txt'

I_hjorth, f_hjorth = unpack(fI_hjorth_file)
I_allen, f_allen = unpack(fI_allen_file)


fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,4),dpi=300)
ax1.plot(t_hjorth,V_hjorth,label='Hjorth et. al.')
ax1.plot(t_allen,V_allen,label='Allen Institute')
ax1.set_xlabel(r'$t$ (ms)', fontsize=12)
ax1.set_ylabel(r'$V$ (mV)', fontsize=12)
ax1.set_title(r'A', loc='left', fontsize=14)
ax1.set_title(r'$I$ = 0.01 nA', fontsize=14)
ax1.axis([-5.7,135.8,-95,51.9]) # 125.8
ax1.legend(loc='lower center', ncol=2, fontsize=12)

ax2.plot(I_hjorth,f_hjorth,label='Hjorth et. al.')
ax2.plot(I_allen,f_allen,label='Allen Institute')
ax2.set_xlabel(r'$I$ (nA)', fontsize=12)
ax2.set_ylabel(r'$f$ (Hz)', fontsize=12)
ax2.set_title(r'B', loc='left', fontsize=14)
ax2.set_title(r'$f-I$', fontsize=14)
ax2.legend(loc='upper left', ncol=1, fontsize=12)
plt.tight_layout()
plt.savefig(figname)
plt.show()


