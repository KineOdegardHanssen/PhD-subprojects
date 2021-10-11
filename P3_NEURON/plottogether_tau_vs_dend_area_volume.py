import numpy as np
import matplotlib.pyplot as plt
import math

dtexp = -10

cm_readoff = 1.0

zoomed = False
somasize = 10
denddiams = [0,0.01,0.1,1,2,5,10,20]
Nd = len(denddiams)
dendlens = [0,1,2,5,10,20,50,100]
Nl = len(dendlens)

areas = []
volumes = []

A      = np.pi*somasize*somasize
Ra     = 100
v_init = -70

idur = 1000   # ms
iamp = -0.1   # nA 
currentfolder = 'current_idur%i_iamp'%idur+str(iamp)+'/dtexp%i/' % dtexp

taus = []
filenames = []
filenames_Cm = []
filenames_R  = []
measured_C   = []
measured_R   = []
somaonly_folder = 'Somaonly/pas/Results/IStim/Soma10/'+currentfolder
filename_nodend = somaonly_folder+'somaonly_cms_idur%i_iamp'%idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_tau_highCms.txt'
filename_nodend_Cm  = somaonly_folder+'somaonly_cms_idur%i_iamp'%idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_Cm_highCms.txt'
filename_nodend_R = somaonly_folder+'somaonly_cms_idur%i_iamp'%idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_Rin_highCms.txt'

## Fix read-in of somaonly
areas.append(0)
volumes.append(0)
######################### Les: ###########################

file = open(filename_nodend,'r')
lines = file.readlines()

for line in lines:
    words = line.split()
    if len(words)>1:
        cm = float(words[0])
        if cm==cm_readoff:
            tau = float(words[1])
            taus.append(tau)
file.close()
    
file_Cm = open(filename_nodend_Cm,'r')
lines   = file_Cm.readlines()
for line in lines:
    words = line.split()
    if len(words)>0:
        cm = float(words[0])
        if cm==cm_readoff:
            measured_C.append(float(words[1])*A)
file_Cm.close()
        
file_R = open(filename_nodend_R,'r')
lines = file_R.readlines()
for line in lines:
    words = line.split()
    if len(words)>1:
        cm = float(words[0])
        if cm==cm_readoff:
            R = float(words[1])
            measured_R.append(R)
file_R.close()


##########################################################




## Antagelig oppdatere dette her: # Inkludere cm_readoff
plotfolder = 'Comparemodels/BAS_vs_somaonly_passive/dtexp%i/' % dtexp
# Area:
plotname_area   = plotfolder + 'BAS_vs_onecomp_cm'+str(cm_readoff)+'_idur%i_iamp'%idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_varyarea_tau'
plotname_C_area = plotfolder + 'BAS_vs_onecomp_cm'+str(cm_readoff)+'_idur%i_iamp'%idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_varyarea_tau_measureC.png'
plotname_R_area = plotfolder + 'BAS_vs_onecomp_cm'+str(cm_readoff)+'_idur%i_iamp'%idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_varyarea_tau_Rin.png'
# Volume:
plotname_volume   = plotfolder + 'BAS_vs_onecomp_cm'+str(cm_readoff)+'_idur%i_iamp'%idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_varyvolume_tau'
plotname_C_volume = plotfolder + 'BAS_vs_onecomp_cm'+str(cm_readoff)+'_idur%i_iamp'%idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_varyvolume_tau_measureC.png'
plotname_R_volume = plotfolder + 'BAS_vs_onecomp_cm'+str(cm_readoff)+'_idur%i_iamp'%idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_varyvolume_tau_Rin.png'
if zoomed==True:
    plotname_area = plotname_area+'_zoomed.png'
    plotname_volume = plotname_volume+'_zoomed.png'
else:
    plotname_area = plotname_area+'.png'
    plotname_volume = plotname_volume+'.png'


for i in range(1,Nd):
    denddiam = denddiams[i]
    for j in range(1,Nl):
        dendlen = dendlens[j]
        area = np.pi*denddiam*dendlen
        volume = np.pi*denddiam*denddiam*dendlen/4.
        areas.append(area)
        volumes.append(volume)
        # Set filenames:
        folder = 'Ball-and-stick models/BAS_passive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam)+'/'+currentfolder
        filename = folder + 'baspass_cms_idur%i_iamp' % idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_tau_highCms.txt'
        filename_Cm = outfilename = folder +'baspass_cms_idur%i_iamp' % idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_Cm_highCms.txt'
        filename_R = outfilename = folder +'baspass_cms_idur%i_iamp' % idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_pas_Rin'
        # tau-file:
        file = open(filename,'r')
        lines = file.readlines()
        
        #print('NUMBER OF LINES:', len(lines), filenames[i]) ### FOR DEBUGGING
        for line in lines:
            words = line.split()
            if len(words)>1:
                cm = float(words[0])
                if cm==cm_readoff:
                    tau = float(words[1])
                    taus.append(tau)
        
        file_Cm = open(filename_Cm,'r')
        lines   = file_Cm.readlines()
        for line in lines:
            words = line.split()
            if len(words)>0:
                cm = float(words[0])
                if cm==cm_readoff:
                    measured_C.append(float(words[1])*A)
        file_Cm.close()
        
        file_R = open(filename_R,'r')
        lines = file_R.readlines()
        for line in lines:
            words = line.split()
            if len(words)>1:
                cm = float(words[0])
                if cm==cm_readoff:
                    R = float(words[1])
                    measured_R.append(R)
        file_R.close()
N = len(areas)
taus = np.array(taus)
areas = np.array(areas)
volumes = np.array(volumes)
measured_C = np.array(measured_C)
measured_R = np.array(measured_R)

argsort_area = areas.argsort()
argsort_volume = volumes.argsort()
areas = np.sort(areas)
volumes = np.sort(volumes)
C_areas = np.zeros(N)
R_areas = np.zeros(N)
tau_areas = np.zeros(N)
C_volumes = np.zeros(N)
R_volumes = np.zeros(N)
tau_volumes = np.zeros(N)

print('argsort_area:',argsort_area)

for i in range(N):
    C_areas[i] = measured_C[argsort_area[i]]
    R_areas[i] = measured_R[argsort_area[i]]
    tau_areas[i] = taus[argsort_area[i]]
    C_volumes[i] = measured_C[argsort_volume[i]]
    R_volumes[i] = measured_R[argsort_volume[i]]
    tau_volumes[i] = taus[argsort_volume[i]]

print('C_areas:',C_areas)

print('areas:',areas)
print('volumes:',areas)
print('taus:',taus)
print('tau_areas:',tau_areas)
print('len(taus):',len(taus))
print('len(tau_areas):',len(tau_areas))
print('len(areas):',len(areas))

############################## Plotting vs area ###################################################
plt.figure(figsize=(6,5))
plt.plot(areas,tau_areas)
plt.xlabel(r'Dendrite area (nm$^2$)')
plt.ylabel(r'$\tau$')
plt.title(r'$\tau$ vs dendrite area')
plt.tight_layout()
plt.savefig(plotname_area)
# 
plt.figure(figsize=(6,5))
plt.plot(areas,C_areas)
plt.xlabel(r'Dendrite area (nm$^2$)')
plt.ylabel(r'Measured $C$')
plt.title(r'Measured $C$ vs dendrite area')
plt.tight_layout()
plt.savefig(plotname_C_area)
#
plt.figure(figsize=(6,5))
plt.plot(areas,R_areas)
plt.xlabel(r'Dendrite area (nm$^2$)')
plt.ylabel(r'$R_{in}$ [$\Omega$]')
plt.title(r'$R_{in}$ vs dendrite area')
plt.tight_layout()
plt.savefig(plotname_R_area)

############################## Plotting vs volume ###################################################
plt.figure(figsize=(6,5))
plt.plot(volumes,tau_volumes)
plt.xlabel(r'Dendrite volume (nm$^3$)')
plt.ylabel(r'$\tau$')
plt.title(r'$\tau$ vs dendrite volume')
plt.tight_layout()
plt.savefig(plotname_volume)
# 
plt.figure(figsize=(6,5))
plt.plot(volumes,C_volumes)
plt.xlabel(r'Dendrite volume (nm$^3$)')
plt.ylabel(r'Measured $C$')
plt.title(r'Measured $C$ vs dendrite volume')
plt.tight_layout()
plt.savefig(plotname_C_volume)
#
plt.figure(figsize=(6,5))
plt.plot(volumes,R_volumes)
plt.xlabel(r'Dendrite volume (nm$^3$)')
plt.ylabel(r'$R_{in}$ [$\Omega$]')
plt.title(r'$R_{in}$ vs dendrite volume')
plt.tight_layout()
plt.savefig(plotname_R_volume)

plt.show()
