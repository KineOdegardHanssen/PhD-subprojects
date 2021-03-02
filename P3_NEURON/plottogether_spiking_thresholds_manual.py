import matplotlib.pyplot as plt
import numpy as np

modelID_perisomatic = 488462965
if modelID_perisomatic==488462965:
    periosomatic_Cm  = np.array([0.5,1.0,2.0,3.0,3.31732779736,4.0,5.0])
    perisomatic_thrs = np.array([0.131,0.132,0.134,0.137,0.138,0.140,0.143])

modelID_allactive = 496497595
if modelID_allactive==496497595:
    allactive_Cm  = np.array([0.01,0.1,0.5,1.0,2.0,3.0,4.0])
    allactive_thrs = np.array([0.196,0.200,0.215,0.227,0.241,0.251,0.259])

ballandstick_Cm = np.array([0.5,1.0,2.0,3.0,4.0,5.0])
ballandstick_thrs = np.array([0.046,0.066,0.096,0.124,0.150,0.175])

folder = 'Comparemodels/%i_%i_ball-and-stick/' % (modelID_perisomatic,modelID_allactive)
plotname = folder+'thr_vs_Cm_%i_%i_ball-and-stick.png' % (modelID_perisomatic,modelID_allactive)


plt.figure(figsize=(6,5))
plt.plot(allactive_Cm, allactive_thrs, '-o', label='Model %i, all active' % modelID_allactive)
plt.plot(periosomatic_Cm, perisomatic_thrs, '-o', label='Model %i, perisomatic' % modelID_perisomatic)
plt.plot(ballandstick_Cm, ballandstick_thrs, '-o', label='Ball-and-stick model')
plt.xlabel(r'Cell capacitance [$\mu$ F/cm$^2$]')
plt.ylabel(r'Spiking threshold $I_{thr}$ [nA]')
plt.title(r'Cell capacitance vs spiking threshold')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig(plotname)