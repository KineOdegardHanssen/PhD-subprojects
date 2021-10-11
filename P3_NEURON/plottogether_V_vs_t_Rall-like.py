import numpy as np
import matplotlib.pyplot as plt

somasize = 10

cm1 = 1.5
cm2 = 0.75

dtexp=-8

Ra = 100
iamp = 0.1 # -0.1
idur = 100
gpas = 0.0003
vpas = -65
epas = vpas
v_init = -65
dendlen  = 1000
denddiam = 2

# File names
infolder_BAS = 'Ball-and-stick models/BAS_passive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/' 
currentfolder = 'current_idur'+str(idur)+'_iamp'+str(iamp)+'/dtexp%i/' % dtexp
infolder_BAS = infolder_BAS+currentfolder 
infilename_BAS = infolder_BAS+'baspass_cm'+str(cm1)+'_idur%.1f_iamp'%idur+str(iamp)+'_Ra%i_gpas'%Ra+str(gpas)+'_vpas' +str(vpas)+'_V.txt' 

folder_onecomp = 'Somaonly/pas/Results/IStim/Soma%i/current_idur%i_iamp'%(somasize,idur)+str(iamp)+'/dtexp%i/' % dtexp
infilename_oc1 = folder_onecomp+'somaonly_cm'+str(cm1)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_gpas'+str(gpas)+'_epas'+str(epas)+'_V.txt' 
infilename_oc2 = folder_onecomp+'somaonly_cm'+str(cm2)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_gpas'+str(gpas)+'_epas'+str(epas)+'_V.txt' 

### Read files: # And normalize?
## Could have looped...
# Read BAS file:
infile_BAS = open(infilename_BAS,'r')
lines = infile_BAS.readlines()

t_BAS = []
V_BAS = []

for line in lines:
    words = line.split()
    if len(words)>1:
        t_BAS.append(float(words[0]))
        V_BAS.append(float(words[1]))

t_BAS = np.array(t_BAS)
V_BAS = np.array(V_BAS)

Vshift_BAS = V_BAS - min(V_BAS)          # Maybe not neccessary in this case
Vnorm_BAS  = Vshift_BAS/max(Vshift_BAS)

infile_BAS.close()

# Read onecomp. file 1:
infile_oc1 = open(infilename_oc1,'r')
lines = infile_oc1.readlines()

t_oc1 = []
V_oc1 = []

for line in lines:
    words = line.split()
    if len(words)>1:
        t_oc1.append(float(words[0]))
        V_oc1.append(float(words[1]))

t_oc1 = np.array(t_oc1)
V_oc1 = np.array(V_oc1)

Vshift_oc1 = V_oc1 - min(V_oc1)
Vnorm_oc1  = Vshift_oc1/max(Vshift_oc1)

infile_oc1.close()


# Read onecomp. file 1:
infile_oc2 = open(infilename_oc2,'r')
lines = infile_oc2.readlines()

t_oc2 = []
V_oc2 = []

for line in lines:
    words = line.split()
    if len(words)>1:
        t_oc2.append(float(words[0]))
        V_oc2.append(float(words[1]))

t_oc2 = np.array(t_oc2)
V_oc2 = np.array(V_oc2)

Vshift_oc2 = V_oc2 - min(V_oc2)
Vnorm_oc2  = Vshift_oc2/max(Vshift_oc2)

infile_oc2.close()

plt.figure(figsize=(6,5))
plt.plot(t_BAS,V_BAS,label=r'Two comp., $\tau=5$ ms')
plt.plot(t_oc1,V_oc1,label=r'One comp., $\tau=5$ ms')
plt.plot(t_oc2,V_oc2,label=r'One comp., $\tau=2.5$ ms')
plt.xlabel('t [ms]')
plt.xlabel('V [mV]')
plt.legend(loc='lower right')


plt.figure(figsize=(6,5))
plt.plot(t_BAS,Vnorm_BAS,label=r'Two comp., $\tau=5$ ms')
plt.plot(t_oc1,Vnorm_oc1,label=r'One comp., $\tau=5$ ms')
plt.plot(t_oc2,Vnorm_oc2,label=r'One comp., $\tau=2.5$ ms')
plt.xlabel('t [ms]')
plt.xlabel('V/max(V)')
plt.legend(loc='lower right')
plt.axis([8,46,-0.01,1.03])

plt.show()

