import numpy
import matplotlib.pyplot as plt

def avg_and_rms(x):
    N = len(x)
    avgx = numpy.mean(x)
    rmsx = 0
    for i in range(N):
        rmsx += (avgx-x[i])**2
    rmsx = numpy.sqrt(rmsx/(N-1))
    return avgx,rmsx

def manual(filename):
    """Manual approach"""

    # Use numpy to read the trace data from the txt file
    data = numpy.loadtxt(filename)
    
    # Time is the first column
    time = data[:, 0]
    # Voltage is the second column
    voltage = data[:, 1]
    
    return time, voltage 

if __name__ == '__main__':
    testmodel  = 488462965 # 478513437 # 478513407 #   ##  ## 485694403 # 489931686 # 480633479 # ## 
    cm         = 1.0
    spikedurat = -40
    idur       = 2000 #100 # ms
    idelay     = 100
    v_init     = -83.7 # mV
    Ra         = 100
    somasize   = 10 # 15 # 
    dendlen    = 1000
    denddiam   = 1
    nsegments  = 200 
    
    varycao = [0.02,200.0]
    
    if testmodel==496497595:
        cm_soma = 1.14805
        cm_dend = 9.98231
        cm_axon = 3.00603
    elif testmodel==497233271:
        cm_soma = 0.783229
        cm_dend = 1.94512
        cm_axon = 8.25387
    elif testmodel==497230463:
        cm_soma = 1.23729
        cm_dend = 2.57923
        cm_axon = 5.02697
    elif testmodel==497233075:
        cm_soma = 1.64168
        cm_dend = 2.83035
        cm_axon = 9.98442
    elif testmodel==488462965:
        cm_soma = 3.31732779736 # Strange values, right?
        cm_dend = 3.31732779736
        cm_axon = 3.31732779736
    elif testmodel==478513407:
        cm_soma = 1.0
        cm_dend = 1.0
        cm_axon = 1.0
    elif testmodel==480633479:
        cm_soma = 0.704866 # 0.704866118957
        cm_dend = 0.704866 # 0.704866118957
        cm_axon = 0.704866 # 0.704866118957
    elif testmodel==478513437:
        cm_soma = 2.34539964752
        cm_dend = 2.34539964752
        cm_axon = 2.34539964752
    elif testmodel==489931686:
        cm_soma = 1.66244903951
        cm_dend = 1.66244903951
        cm_axon = 1.66244903951
    elif testmodel==485694403:
        cm_soma = 0.763348896
        cm_dend = 0.763348896
        cm_axon = 0.763348896
    
    if testmodel==480633479:
        v_init = -96.8#-83.7#-90#-86.5# # Have a test here too
    elif testmodel==496497595:
        v_init = -86.5
    elif testmodel==488462965:
        v_init = -86.5 # Maybe I should have changed this...
    elif testmodel==497230463:
        v_init = -90
    elif testmodel==497233075:
        v_init = -90
    elif testmodel==478513437:
        v_init = -86.8
    elif testmodel==478513407:
        v_init = -83.7
    elif testmodel==497233271:
        v_init = -90
    elif testmodel==489931686:
        v_init = -95.7
    elif testmodel==485694403:
        v_init = -88.8
    ###########################################################
    
    iamp  = 0.3
    
    Vs = []
    plt.figure(figsize=(6,5))
    # Set names
    outfolder = 'figures/%i/' % testmodel
    for cao in varycao:
        namestring = '_cao'+str(cao)
        infolder = 'figures/%i/current_idur%i_iamp' % (testmodel,idur) + str(iamp)+'/'
        filename = infolder+namestring+"_changecmf" + str(cm) + "_everywhere_vinit"+str(v_init)+"_addedRa.txt"
        t, V = manual(filename)
        plt.plot(t,V)
        Vs.append(V)
    plt.xlabel('t (ms)')
    plt.ylabel('V (mV)')
    plt.show()
    
    print('Vs[1]-Vs[0]:',Vs[1]-Vs[0])
    print('sum(Vs[1]-Vs[0]):',sum(Vs[1]-Vs[0]))
    print('max(Vs[1]-Vs[0]):',max(Vs[1]-Vs[0]))
    print('min(Vs[1]-Vs[0]):',min(Vs[1]-Vs[0]))
        
    