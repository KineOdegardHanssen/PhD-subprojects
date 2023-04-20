import numpy as np
from matplotlib import gridspec
import matplotlib.pyplot as plt

def unpack(filename):
    data = np.loadtxt(filename)

    # Time is the first column
    x = data[:, 0]
    # Voltage is the second column
    y = data[:, 1]
    return x, y

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

figname = 'choosemodels_trace_fI_Cain.png'

# Trace file
folderstart_Hjorth = 'Hjorth2020/FS_Na_K/'
allenfile = 'Allen_NaV_Kv2like/Results/Soma10/current_idur100_iamp0.01/Epasshift-20_gNaV0.5p_gKv2like0.5p_somaonly_cm1.0_idur100_iamp0.01_dtexp-7_V.txt'
hjorthfile = folderstart_Hjorth+'Results/Soma10/current_idur100_iamp0.01/somaonly_cm1.0_idur100_iamp0.01_dtexp-7_vinit-86_trec-600.0_V.txt'

t_allen, V_allen = unpack(allenfile)
t_hjorth, V_hjorth = unpack(hjorthfile)

# fI-curve
fI_hjorth_file = folderstart_Hjorth+'/Results/Soma10/somaHjorth_idur1000_varyiamp_manual_cm1.0_Nspikes_vs_I_s500.txt'
fI_allen_file  = 'Allen_NaV_Kv2like/Results/Soma10/somaAllen_idur1000_varyiamp_manual_cm1.0_Epasshift-20_gNaV0.5p_gKv2like0.5p_Nspikes_vs_I_s500.txt'

I_hjorth, f_hjorth = unpack(fI_hjorth_file)
I_allen, f_allen = unpack(fI_allen_file)

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

folderstart_Hjorth = folderstart_Hjorth + 'CaHVA_Allen/'
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
    filename = folderstart_Hjorth+folderbase+'_gCaHVA'+str(gcahva_this)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V_cainfo.txt'
    
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

gsks = [1]

filename_sk1  = folderstart_Hjorth+'SK_Allen/'+folderbase+'_gSK'+str(gsk)+'p_gCaHVA'+str(gcahva)+'p_somaonly_cm1.0_idur'+str(idur)+'_iamp'+str(iamp)+'_dtexp-7_vinit-86_trec-600.0_V_eca.txt'
print('filename:',filename_sk1)


t_sk1, V_sk1, eca_sk1, cai_sk1, cao_sk1, I_SK_sk1, I_Ca_HVA_sk1, g_SK_sk1, g_Ca_HVA_sk1 = unpack_sk(filename_sk1)

g_SK_sk1 = g_SK_sk1/float(gsk*gsk_base)

gcahva_base = 0.00014931667074610222
g_Ca_HVA_sk1  = g_Ca_HVA_sk1/float(gcahva*gcahva_base)


fig = plt.figure(figsize=(10,6),dpi=300)
    
gs = gridspec.GridSpec(2, 4)
    
ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[1, 0:2])
ax4 = plt.subplot(gs[1, 2:4])

ax1.set_title('A',loc='left')
ax2.set_title('B',loc='left')
ax3.set_title('C',loc='left')
ax4.set_title('D',loc='left')

ax1.plot(t_hjorth,V_hjorth)
ax1.plot(t_allen,V_allen)
ax1.set_xlabel(r'$t$ (ms)', fontsize=12)
ax1.set_ylabel(r'$V$ (mV)', fontsize=12)
ax1.set_title(r'A', loc='left', fontsize=14)
ax1.set_title(r'$I$ = 0.01 nA', fontsize=14)
#ax1.axis([-5.7,135.8,-95,51.9]) # 125.8
#ax1.legend(loc='lower center', ncol=2, fontsize=12)

ax2.plot(I_hjorth,f_hjorth,label='Hjorth et. al.')
ax2.plot(I_allen,f_allen,label='Allen Institute')
ax2.set_xlabel(r'$I$ (nA)', fontsize=12)
ax2.set_ylabel(r'$f$ (Hz)', fontsize=12)
ax2.set_title(r'B', loc='left', fontsize=14)
ax2.set_title(r'$f-I$', fontsize=14)
ax2.legend(loc='upper left', ncol=1, fontsize=12)

ax3.plot(t1, cai1)
ax3.set_xlabel(r'$t$ (ms)', fontsize=12)
ax3.set_title(r'C', loc='left', fontsize=14)
ax3.set_ylabel(r'$[\mathregular{Ca}^{2+}]_\mathregular{in}$ (mM)', fontsize=12)
ax3.set_title(r'$[\mathregular{Ca}^{2+}]_\mathregular{in}$, CaHVA, $I$ = 0.02 nA', fontsize=14)

ax4.plot(t_sk1, cai_sk1)
ax4.set_xlabel(r'$t$ (ms)', fontsize=12)
ax4.set_title(r'D', loc='left', fontsize=14)
ax4.set_ylabel(r'$[\mathregular{Ca}^{2+}]_\mathregular{in}$ (mM)', fontsize=12)
ax4.set_title(r'$[\mathregular{Ca}^{2+}]_\mathregular{in}$, CaHVA and SK, $I$ = 0.02 nA', fontsize=13.5)

plt.tight_layout()

plt.savefig(figname)
plt.show()


