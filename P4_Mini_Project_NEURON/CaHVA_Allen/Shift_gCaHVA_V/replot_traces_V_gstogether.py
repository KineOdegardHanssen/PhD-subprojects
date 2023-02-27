import numpy
import matplotlib.pyplot as plt
from matplotlib import gridspec

spikedurat = -40
idur       = 1000 #100 # ms
idelay     = 10
v_init     = -86 # mV
Ra         = 100
somasize   = 10 # 15 # 
dtexp      = -7
t_before_rec = -600.

vary_naf     = False
vary_kaf     = False
vary_Ca_HVA  = True # False
vary_gpas    = False # 

cm       = 1.0
iamp     = 0.005
gcahvas  = [0.1,1.0]#[0.1,0.5,1.0,1.5]
vshifts  = [-25,-5,0,5]#[-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50]
Ng       = len(gcahvas)
Nv       = len(vshifts)

colors = []
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))

outfolder    = 'Compare/'
infolder_end = 'Results/Soma%i/' % somasize+'current_idur%i_iamp'%idur+str(iamp)+'/'


plotname = outfolder + 'selectedvshift_selectedcahva_iamp'+str(iamp)+'_subplotV_gt_traces_Ng%i.png' % Ng
fig = plt.figure(figsize=(18,18),dpi=300)

gs = gridspec.GridSpec(2, 4)

ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[1, 0:2])
ax4 = plt.subplot(gs[1, 2:4])

#fig.suptitle(r'Properties',fontsize=20)

ax1.set_title(r'A',loc='left',fontsize=18)
ax2.set_title(r'B',loc='left',fontsize=18)
ax3.set_title(r'C',loc='left',fontsize=18)
ax4.set_title(r'D',loc='left',fontsize=18)
ax1.set_title(r'$V_{\mathregular{Shift}}$=%i' % vshifts[0],fontsize=18)
ax2.set_title(r'$V_{\mathregular{Shift}}$=%i' % vshifts[1],fontsize=18)
ax3.set_title(r'$V_{\mathregular{Shift}}$=%i' % vshifts[2],fontsize=18)
ax4.set_title(r'$V_{\mathregular{Shift}}$=%i' % vshifts[3],fontsize=18)

ax1.set_xlabel(r'$t$ [ms]')    
ax1.set_ylabel(r'$V$ [mV]')
ax2.set_xlabel(r'$t$ [ms]')
ax2.set_ylabel(r'$V$ [mV]')
ax3.set_xlabel(r'$t$ [ms]')
ax3.set_ylabel(r'$V$ [mV]')
ax4.set_xlabel(r'$t$ [ms]')
ax4.set_ylabel(r'$V$ [mV]')

times    = []
voltages = []
for i in range(Nv):
    vshift = vshifts[i]
    infolder = 'Vshift_'+str(vshift)+'/'+infolder_end
    
    times_this    = []
    voltages_this = []
    
    for j in range(Ng):
        gcahva = gcahvas[j]
        print('Step ', j+1, ' of', Ng, 'cm:', cm, 'gcahva:', gcahva)    
        namestring = ''
        thelabel = ''
        if vary_naf==True:
            namestring = namestring + '_gnaf'+str(gnaf)+'p'
            #thelabel   = thelabel + r'$\bar{g}_\mathregular{naf}$=%.1f, ' % gnaf
        if vary_kaf==True:
            namestring = namestring + '_gkaf'+str(gkaf)+'p'
            #thelabel   = thelabel + r'$\bar{g}_\mathregular{kaf}$=%.1f, ' % gkaf
        if vary_Ca_HVA==True:
            namestring = namestring + '_gCaHVA'+str(gcahva)+'p'
            thelabel   = thelabel + r'$\bar{g}_\mathregular{CaHVA}$=%.1f' % gcahva
        if vary_gpas==True: 
            namestring = namestring + '_gpas'+str(gpas)+'p'
            #thelabel   = thelabel + r'$\bar{g}_\mathregular{L}$=%.1f, ' % gpas
        namestring = namestring +'_'
        filename = infolder+namestring+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_V.txt'
        
        # Extract the data
        data = numpy.loadtxt(filename)# Time is the first column
        time = data[:, 0]
        # Voltage is the second column
        voltage = data[:, 1]
        
        # Append data
        times_this.append(time)
        voltages_this.append(voltage)
    times.append(times_this)
    voltages.append(voltages_this)


times_this    = times[0]
voltages_this = voltages[0]
for i in range(Ng):
    time    = times_this[i]
    voltage = voltages_this[i]
    ax1.plot(time, voltage, label=r'%s$\bar{g}_{\mathregular{CaHVA}}$' % str(gcahvas[i]),linewidth=2,color=colors[i])


times_this    = times[1]
voltages_this = voltages[1]
for i in range(Ng):
    time    = times_this[i]
    voltage = voltages_this[i]
    ax2.plot(time, voltage, label=r'%s$\bar{g}_{\mathregular{CaHVA}}$' % str(gcahvas[i]),linewidth=2,color=colors[i])


times_this    = times[2]
voltages_this = voltages[2]
for i in range(Ng):
    time    = times_this[i]
    voltage = voltages_this[i]
    ax3.plot(time, voltage, label=r'%s$\bar{g}_{\mathregular{CaHVA}}$' % str(gcahvas[i]),linewidth=2,color=colors[i])


times_this    = times[3]
voltages_this = voltages[3]
for i in range(Ng):
    time    = times_this[i]
    voltage = voltages_this[i]
    ax4.plot(time, voltage, label=r'%s$\bar{g}_{\mathregular{CaHVA}}$' % str(gcahvas[i]),linewidth=2,color=colors[i])

ax1.legend(loc='upper left')#,ncol=1)
ax2.legend(loc='upper left')#,ncol=1)
ax3.legend(loc='upper left')#,ncol=1)
ax4.legend(loc='upper left')#,ncol=1)
plt.tight_layout()#(rect=[0, 0.03, 1, 0.95])
plt.savefig(plotname)
    
