import os
from os.path import join
import sys
import matplotlib.pyplot as plt
import json
import neuron
import time as tm
import numpy as np
from matplotlib import gridspec

iamp = 0.006
idur = 1000
dtexp = -7
v_init = -86
somasize = 10
cm_factor = 1.0
t_before_rec = -600.0
conc_at_halfopen = 0.00043287612810830614

gcahva = 0.2
gsk    = 1.0

namestring = ''
namestring = namestring + '_gSK'+str(gsk)+'p'
namestring = namestring + '_gCaHVA'+str(gcahva)+'p'
namestring = namestring +'_'

folder = 'Results/Soma%i/current_idur%i_iamp'%(somasize,idur)+str(iamp)+'/'
if os.path.exists(folder)==False:
    os.mkdir(folder)
#t, v, eca, cai, cao
filename = folder+namestring+'somaonly_cm'+str(cm_factor)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_V_eca.txt'

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


#fig, (ax1, ax2, ax3) = plt.subplots(3,1,figsize=(5,11),dpi=300)
fig = plt.figure(figsize=(10,8),dpi=300)#(figsize=(8,3),dpi=300)
    
gs = gridspec.GridSpec(4, 4)

ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[1, 0:2])
ax3 = plt.subplot(gs[2, 0:2])
ax4 = plt.subplot(gs[3, 0:2])
ax5 = plt.subplot(gs[0, 2:4])
ax6 = plt.subplot(gs[1, 2:4])
ax7 = plt.subplot(gs[2, 2:4])
ax8 = plt.subplot(gs[3, 2:4])


ax1.plot(t,v)
ax1.axvline(x=100,color='k',linestyle=':',linewidth=0.75)
ax1.axvline(x=1100,color='k',linestyle=':',linewidth=0.75)
#ax1.set_xlabel('V (mV)')
ax1.set_ylabel(r'$V$ (mV)',fontsize=12)
ax1.set_title(r'$I=$ %s nA' % str(iamp),fontsize=16)
ax1.set_title('A', loc='left',fontsize=18)

ax2.plot(t,eca,color='k')
ax2.axvline(x=100,color='k',linestyle=':',linewidth=0.75)
ax2.axvline(x=1100,color='k',linestyle=':',linewidth=0.75)
#ax1.set_xlabel('V (mV)')
ax2.set_ylabel(r'$E_\mathregular{Ca}$',fontsize=12)
ax2.set_title(r'$E_\mathregular{Ca}$',fontsize=16)
ax2.set_title('A', loc='left',fontsize=18)

ax3.plot(t,cai,color='tab:brown')
ax3.axvline(x=100,color='k',linestyle=':',linewidth=0.75)
ax3.axhline(y=conc_at_halfopen,color='k',linestyle='--',linewidth=0.75)
ax3.axvline(x=1100,color='k',linestyle=':',linewidth=0.75)
#ax2.set_xlabel('t (ms)',fontsize=12)
ax3.set_ylabel(r'Concentration (mM)',fontsize=12)
ax3.set_title(r'$\left[\mathregular{Ca}^{2+}\right]_\mathregular{in}$',fontsize=16)
ax3.set_title('B', loc='left',fontsize=18)

ax4.plot(t,cao,color='tab:brown')
ax4.axvline(x=100,color='k',linestyle=':',linewidth=0.75)
ax4.axvline(x=1100,color='k',linestyle=':',linewidth=0.75)
ax4.set_xlabel(r'$t$ (ms)',fontsize=12)
ax4.set_ylabel(r'Concentration (mM)',fontsize=12)
ax4.set_title(r'$\left[\mathregular{Ca}^{2+}\right]_\mathregular{out}$',fontsize=16)
ax4.set_title('C', loc='left',fontsize=18)

ax5.plot(t,I_SK,color='tab:gray')
#ax1.set_xlabel('V (mV)')
ax5.axvline(x=100,color='k',linestyle=':',linewidth=0.75)
ax5.axvline(x=1100,color='k',linestyle=':',linewidth=0.75)
ax5.set_ylabel(r'$I_\mathregular{SK}$ (nA)',fontsize=12)
ax5.set_title(r'$I_\mathregular{SK}$',fontsize=16)
ax5.set_title('D', loc='left',fontsize=18)

ax6.plot(t,I_Ca_HVA,color='tab:gray')
#ax1.set_xlabel('V (mV)')
ax6.axvline(x=100,color='k',linestyle=':',linewidth=0.75)
ax6.axvline(x=1100,color='k',linestyle=':',linewidth=0.75)
ax6.set_ylabel(r'$I_\mathregular{CaHVA}$ (nA)',fontsize=12)
ax6.set_title(r'$I_\mathregular{CaHVA}$',fontsize=16)
ax6.set_title('E', loc='left',fontsize=18)

ax7.plot(t,g_SK,color='tab:purple')
#ax1.set_xlabel('V (mV)')
ax7.axvline(x=100,color='k',linestyle=':',linewidth=0.75)
ax7.axvline(x=1100,color='k',linestyle=':',linewidth=0.75)
ax7.set_ylabel(r'$g_\mathregular{SK}$ (S/cm$^2$)',fontsize=12)
ax7.set_title(r'$g_\mathregular{SK}$',fontsize=16)
ax7.set_title('F', loc='left',fontsize=18)

ax8.plot(t,g_Ca_HVA,color='tab:purple')
#ax1.set_xlabel('V (mV)')
ax8.axvline(x=100,color='k',linestyle=':',linewidth=0.75)
ax8.axvline(x=1100,color='k',linestyle=':',linewidth=0.75)
ax8.set_xlabel(r'$t$ (ms)',fontsize=12)
ax8.set_ylabel(r'$g_\mathregular{CaHVA}$ (S/cm$^2$)',fontsize=12)
ax8.set_title(r'$g_\mathregular{CaHVA}$',fontsize=16)
ax8.set_title('G', loc='left',fontsize=18)

plt.tight_layout()
plt.savefig('Results/Soma%i/Ca-properties/Ca_E_and_concentrations_iamp'%somasize+str(iamp)+'_idur'+str(idur)+'.png')
plt.show()


