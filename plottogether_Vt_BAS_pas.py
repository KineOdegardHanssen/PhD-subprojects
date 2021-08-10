import numpy as np
import matplotlib.pyplot as plt
import math

padding = True
# loop de loop? # Yes.
somasize = 10
dendlen  = 1
denddiam = 20
L      = somasize
diam   = somasize
A      = np.pi*L*diam
Ra     = 100
v_init = -70
long   = True
    
# Values of membrane capacitance (parameters):
cms = [0.01,0.1,0.5,1.0,1.5,2.0,5.0,15.0]
Ncm = len(cms)
# To choose from: [0.01,0.1,0.5,0.8,1.0,1.2,1.5,2.0,3.0,5.0,7.0,10.0,15.0]

cmsstring = str(cms[0])
for i in range(1,Ncm):
    cmsstring = cmsstring + '_'+str(cms[i]) # Too long to use

# Change current:
idur = 1000   # ms
iamp = -0.1   # nA 
idelay = 100 

print('Ncm:',Ncm)

folder = 'Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam) +'/current_idur%i_iamp'% idur+str(iamp)+'/'
plotname = folder +'baspass_cms_idur%i_iamp' % idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_Vt.png'

Vlist = []

Ncm = 0
for cm in cms:
    infilename  = folder+'baspass_cm'+str(cm)+'_idur%.1f_iamp' % idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_V.txt'
    infile = open(infilename,'r')
    lines = infile.readlines()
    N = len(lines) # Do I need this?
    
    time = []
    V    = []
    
    for line in lines:
        words = line.split()
        time.append(float(words[0]))
        V.append(float(words[1]))
    infile.close()
    
    Vlist.append(V)
    Vmax = max(V)
    Vmin = min(V)
    Ncm+=1
    

print('Ncm:',Ncm)

plt.figure(figsize=(6,5))
for i in range(Ncm):
    plt.plot(time,Vlist[i],label=r'$C_m=$%s' % str(cms[i]))
plt.xlabel(r'$t$ [ms]')
plt.ylabel(r'$V$ [mV]')
plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f' % (idur,iamp))
plt.legend(loc='upper right')
plt.axis([idelay-5,200,Vmin-padding,Vmax+padding])
plt.tight_layout()
plt.savefig(plotname)
plt.show()