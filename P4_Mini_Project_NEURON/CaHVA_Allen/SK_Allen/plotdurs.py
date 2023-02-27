import numpy
import matplotlib.pyplot as plt

spikedurat = -40
idur       = 1000 #100 # ms
idelay     = 10
v_init     = -86.8 # mV
Ra         = 100
somasize   = 10 # 15 # 
dtexp      = -7
    
outfolder = 'Results/Soma%i/' % somasize
    
cm = 1.0
iamps =  [0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03,0.031,0.032,0.033,0.034,0.035,0.036,0.037,0.038,0.039,0.04,0.041,0.042,0.043,0.044,0.045,0.046,0.047,0.048,0.049,0.05,0.051,0.052,0.053,0.054,0.055,0.056,0.057,0.058,0.059,0.06]

Namps = len(iamps)

for j in range(Namps):
    durs  = []
    times = []
    iamp = iamps[j]
    print('Step ', j+1, ' of', Namps, 'cm:', cm, 'iamp:', iamp)
    infolder = outfolder + 'current_idur%i_iamp'%idur+str(iamp)+'/'
    figname = outfolder+'Spikedurs/somaHjorth_idur%i_iamp'% (idur)+str(iamp)+'.png'
    ### Spike duration ###
    filename = infolder+'somaHjorth_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmf'+str(cm)+'_alldurs.txt'
    try:
        infile = open(filename,'r')
    except:
        continue
    line = infile.readline()
    words = line.split()
    for word in words:
        durs.append(float(word))
    infile.close()
    ### Spike times ###
    tfilename = infolder+'somaHjorth_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmf'+str(cm)+'_peaktimes.txt'
    tfile = open(tfilename,'r')
    line = tfile.readline()
    words = line.split()
    for word in words:
        times.append(float(word))
    tfile.close()
    ###
    times = times[1:-1]
    print('len(times):',len(times))
    print('len(durs):',len(durs))
    if len(times)==len(durs):
        plt.figure(figsize=(6,5))
        plt.plot(times,durs)
        plt.xlabel('Time [ms]')
        plt.ylabel('Voltage [mV]')
        plt.title('I=%s nA' % str(iamp))
        plt.tight_layout()
        plt.savefig(figname)
