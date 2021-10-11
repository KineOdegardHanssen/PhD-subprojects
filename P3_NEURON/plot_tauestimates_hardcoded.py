import matplotlib.pyplot as plt
import numpy as np

cms = np.array([0.5,1.0,1.5,2.0,3.0,5.0])
efeltau = np.array([2.863,4.925,8.085,11.943,19.180,12.603])
readtau = np.array([0.8777,1.63,2.41,3.22,4.937,8.658])

plt.figure(figsize=(6,5))
plt.plot(cms,efeltau, label='eFEL')
plt.plot(cms,readtau, label='Manually')
plt.xlabel(r'$C_m$ soma parameter [$\mu$F/cm$^2$]')
plt.ylabel(r'$\tau$ [ms]')
plt.title(r'Different estimates of $\tau$')
plt.legend(loc='upper left')
plt.savefig('tauestimates.png')

plt.plot(cms,readtau, 'o',label='Manually')
plt.show()