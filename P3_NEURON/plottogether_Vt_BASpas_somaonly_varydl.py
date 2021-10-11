import numpy as np
import matplotlib.pyplot as plt
import math

dtexp = -8

plotvsC  = False # True # 
plotvsd  = True # False # 
plotvsl  = False # True # 
paddit   = True
zerobool = False # DO NOT TOUCH!
if paddit==True:
    padding = 2
else:
    padding = 0
# loop de loop? # Yes.
somasize = 10
dendlen  = 1000
denddiam = 1
cm     = 1.0
cmb    = 1.0
Ra     = 100
v_init = -65
gpas   = 0.0003
vpas   = -65
long   = True
soma_v_init = -65
dendlen  = 1000
denddiam = 2

time0 = []
    
# Values of membrane capacitance (parameters):
if plotvsC==True:
    cms = [0.01,0.1,0.5,1.0,1.5,2.0,5.0,15.0]
    Npl = len(cms)
    # To choose from: [0.01,0.1,0.5,0.8,1.0,1.2,1.5,2.0,3.0,5.0,7.0,10.0,15.0]
elif plotvsd==True: # Include somaonly in these? A lot of work
    denddiams = [0,0.01,0.1,1,2,5,10,20]
    # To choose from: denddiams = [0,0.01,0.1,1,2,5,10,20]
    Npl = len(denddiams)
elif plotvsl==True:
    dendlens = [0,1,2,5,10,20,50,100]
    # To choose from: dendlens = [1,2,5,10,20,50,100]
    Npl = len(dendlens)

#cmsstring = str(cms[0])
#for i in range(1,Ncm):
#    cmsstring = cmsstring + '_'+str(cms[i]) # Too long to use

# Change current:
idur = 100   # ms
iamp = -0.1   # nA 
idelay = 10

taus = []
tauheights = []
tauheights_norm = []

print('Npl:',Npl)

### Setting file name
if plotvsC==True:
    folder = 'Ball-and-stick models/BAS_passive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam) +'/current_idur%i_iamp'% idur+str(iamp)+'/'
    plotname = folder +'baspass_cms_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+str(v_init)+'_pas_Vt.png'
    plotname_2 = folder +'baspass_cms_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+'_gpas'+str(gpas)+'_vpas'+str(vpas)+'_vinit'+str(v_init)+'_pas_Vt_zoomed.png'
    plotname_norm = folder +'baspass_cms_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+'_gpas'+str(gpas)+'_vpas'+str(vpas)+'_vinit'+str(v_init)+'_pas_Vt_normalized.png'
    plotname_crosses = folder +'baspass_cms_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+'_gpas'+str(gpas)+'_vpas'+str(vpas)+'_vinit'+str(v_init)+'_pas_Vt_withtaus.png'
    plotname_withfit = folder +'baspass_cms_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+'_gpas'+str(gpas)+'_vpas'+str(vpas)+'_vinit'+str(v_init)+'_pas_Vt_withfit.png'
    Vcut = 200
elif plotvsd==True:
    folder = 'Comparemodels/BAS_vs_somaonly_passive/dtexp%i/Varyds/' % dtexp
    plotname = folder +'baspass_cm'+str(cm)+'_cmb'+str(cmb)+'_varyds_dendlen'+str(dendlen)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+str(v_init)+'_pas_Vt.png'
    plotname_2 = folder +'baspass_cm'+str(cm)+'_cmb'+str(cmb)+'_varyds_dendlen'+str(dendlen)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+'_gpas'+str(gpas)+'_vpas'+str(vpas)+'_vinit'+str(v_init)+'_pas_Vt_zoomed.png'
    plotname_norm = folder +'baspass_cm'+str(cm)+'_cmb'+str(cmb)+'_varyds_dendlen'+str(dendlen)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+'_gpas'+str(gpas)+'_vpas'+str(vpas)+'_vinit'+str(v_init)+'_pas_Vt_normalized.png'
    plotname_crosses = folder +'baspass_cm'+str(cm)+'_cmb'+str(cmb)+'_varyds_dendlen'+str(dendlen)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+'_gpas'+str(gpas)+'_vpas'+str(vpas)+'_vinit'+str(v_init)+'_pas_Vt_withtaus.png'
    plotname_withfit = folder +'baspass_cm'+str(cm)+'_cmb'+str(cmb)+'_varyds_dendlen'+str(dendlen)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+'_gpas'+str(gpas)+'_vpas'+str(vpas)+'_vinit'+str(v_init)+'_pas_Vt_withfit.png'
    Vcut = idelay+25 # Suitable for cm=1.0. # COULD NEED TO UPDATE THIS. Possibly test for the value of cm.
    Vcut_zoom = idelay+5 # Also need to beware of this
elif plotvsl==True:
    folder = 'Comparemodels/BAS_vs_somaonly_passive/dtexp%i/Varyls/' % dtexp
    plotname = folder +'bashh_cm'+str(cm)+'_varyls_denddiam'+str(denddiam)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+'_gpas'+str(gpas)+'_vpas'+str(vpas)+'_pas_Vt.png'
    plotname_2 = folder +'baspass_cm'+str(cm)+'_cmb'+str(cmb)+'_varyls_denddiam'+str(denddiam)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+'_gpas'+str(gpas)+'_vpas'+str(vpas)+'_pas_Vt_zoomed.png'
    plotname_norm = folder +'baspass_cm'+str(cm)+'_cmb'+str(cmb)+'_varyls_denddiam'+str(denddiam)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+'_gpas'+str(gpas)+'_vpas'+str(vpas)+'_pas_Vt_normalized.png'
    plotname_crosses = folder +'baspass_cm'+str(cm)+'_cmb'+str(cmb)+'_varyls_denddiam'+str(denddiam)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+'_gpas'+str(gpas)+'_vpas'+str(vpas)+'_pas_Vt_withtaus.png'
    plotname_withfit = folder +'baspass_cm'+str(cm)+'_cmb'+str(cmb)+'_varyls_denddiam'+str(denddiam)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i'%Ra+'_gpas'+str(gpas)+'_vpas'+str(vpas)+'_pas_Vt_withfit.png'
    Vcut = idelay+25 # Suitable for cm=1.0. # COULD NEED TO UPDATE THIS. Possibly test for the value of cm.
    Vcut_zoom = idelay+5 # Also need to beware of this

Vlist   = []
fitlist = []
Vnormlist  = []
Vshiftlist = []
Vnormmin   = 1

for i in range(Npl):
    ### Setting file name according to what you want to vary
    if plotvsC==True:
        cm = cms[i]
        infilename  = folder+'baspass_cmf'+str(cmfac)+'_idur%i_iamp' % idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_V.txt'
    elif plotvsd==True:
        denddiam = denddiams[i]
        folder = 'Ball-and-stick models/BAS_passive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam) +'/current_idur%i_iamp'% idur+str(iamp)+'/dtexp%i/' % dtexp
        infilename  = folder+'baspass_cm'+str(cm)+'_cmb'+str(cmb)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i'%Ra+'_gpas'+str(gpas)+'_vpas'+str(vpas)+'_sprx_V.txt' 
        print('dendlen:', dendlen, ', denddiam:',denddiam)
        # Testing for denddiam==0. If so, get soma data
        if denddiam==0:
            somaonly_folder = 'Somaonly/pas/Results/IStim/Soma10/current_idur%i_iamp'% idur+str(iamp)+'/dtexp%i/' % dtexp
            infilename = somaonly_folder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+'_gpas'+str(gpas)+'_epas'+str(vpas)+'_V.txt' 
            zerobool = True
    elif plotvsl==True:
        dendlen = dendlens[i]
        folder = 'Ball-and-stick models/BAS_passive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam) +'/current_idur%i_iamp'% idur+str(iamp)+'/dtexp%i/' % dtexp
        infilename  = folder+'baspass_cm'+str(cm)+'_cmb'+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_gpas' %Ra+str(g_init)+'_vpas'+str(vpas)+'_sprx_V.txt' 
        print('dendlen:',dendlen, ', denddiam:',denddiam)
        # Testing for dendlen==0. If so, get soma data
        if dendlen==0:
            somaonly_folder = 'Somaonly/pas/Results/IStim/Soma10/current_idur%i_iamp'% idur+str(iamp)+'/dtexp%i/' % dtexp
            infilename = somaonly_folder+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_ena'%Ra+str(ena)+'_ek' +str(ek)+'_el' +str(el)+'_gnabar'+str(gnabar)+'_gkbar'+str(gkbar)+'_gl'+str(gl)+'_Ra%i_vinit'%Ra+str(soma_v_init)+'_gpas'+str(gpas)+'_epas' +str(vpas)+'_V.txt' 
            zerobool = True
    
    ### Opening and reading the file
    infile = open(infilename,'r')
    lines = infile.readlines()
    N = len(lines) # Do I need this?
    
    times = []
    V     = []
    
    for line in lines:
        words = line.split()
        times.append(float(words[0]))
        V.append(float(words[1]))
    infile.close()
    
    Vlist.append(V)
    if plotvsC==True:
        if i==Npl-1:
            Vmax = max(V)
            Vmin = min(V)
    elif plotvsd==True or plotvsl==True:
        if i==0:
            Vmax = max(V)
            Vmin = min(V)
    
    Vmaxf = V[0] #max(V)
    Vminf = min(V)
    dV = Vmaxf-Vminf
    dV1e = dV*math.exp(-1)
    V1e  = Vminf+dV1e
    #print('Vmaxf:',Vmaxf)
    #print('Vminf:',Vminf)
    #print('min(V):',min(V))
    #print('V1e:',V1e)
    #print('dV:',dV)
    #print('dV1e:',dV1e)

    # Will perform linear interpolation near the crossing with V1e.
    tbef = 0
    taft = 0
    Vbef = 0
    Vaft = 0
    for j in range(len(V)):
        if V[j]-V1e<0:
            tbef = times[j-1]
            taft = times[j]
            Vbef = V[j-1]
            Vaft = V[j]
            break
    
    #print('i:',i, '; tbef=',tbef, '; taft=',taft)
    a = (Vaft-Vbef)/(taft-tbef)
    idelay_calc = idelay
    if zerobool==True:
        idelay_calc = 10 # Inconsistency in idelay, somaonly vs. BAS
    t_interpolated = (V1e-Vbef+a*tbef)/a-idelay_calc
    
    tau = t_interpolated
    
    #print('dV+Vmaxf:',dV+Vmaxf)
     # Could do fits later...
    Nf = 100*idur
    timecut = np.linspace(0,idur,Nf)
    fit = np.zeros(Nf)
    for j in range(Nf):
        fit[j] = dV*math.exp(-timecut[j]/tau)+Vminf
    timecut = np.linspace(idelay,idur+idelay,Nf)
    fitlist.append(fit)
    
    #print('i:',i, '; tau:', tau, '; t_interpolated:', t_interpolated)
    taus.append(tau+idelay)
    tauheights.append(V1e) # Should I have normalized first?
    tauheights_norm.append(math.exp(-1))
    
    V = np.array(V)
    Vshift = V-min(V)
    Vnorm  = Vshift/max(Vshift)
    if zerobool==True: # Not sure this is necessary. The bug was somewhere else?
        time0  = np.array(times)
        V0     = np.array(V)
        V0shift = V0-min(V0)
        V0norm = V0shift/max(V0shift)
        Vnorm  = V0norm 
        Vnormlist.append(V0norm)
        Vshiftlist.append(V0shift)
    else:
        Vnormlist.append(Vnorm)
        Vshiftlist.append(Vshift)
    
    zerobool=False
    
    minVnow = min(Vnorm)
    if minVnow<Vnormmin:
        Vnormmin = minVnow

taus = np.array(taus)
tauheights = np.array(tauheights)
tauheights_norm = np.array(tauheights_norm)

print('len(taus):',len(taus))
print('taus:',taus)

plt.figure(figsize=(6,5))
if plotvsC==True:
    plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f, d=%.2f, l=%i' % (idur,iamp,denddiam,dendlen))
    for i in range(Npl):
        plt.plot(times,Vlist[i],label=r'$C_m=$%s' % str(cms[i]))
elif plotvsd==True:
    plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f, Cm=%.2f, l=%i' % (idur,iamp,cm,dendlen))
    plt.plot(time0,V0,label=r'$d=$%s' % str(denddiams[0]))
    for i in range(1,Npl):
        plt.plot(times,Vlist[i],label=r'$d=$%s' % str(denddiams[i]))
elif plotvsl==True:
    plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f, Cm=%.2f, d=%.2f' % (idur,iamp,cm,denddiam))
    plt.plot(time0,V0,label=r'$l=$%s' % str(dendlens[0]))
    for i in range(1,Npl):
        plt.plot(times,Vlist[i],label=r'$l=$%s' % str(dendlens[i]))
plt.xlabel(r'$t$ [ms]')
plt.ylabel(r'$V$ [mV]')
plt.legend(loc='upper right')
plt.axis([idelay-5,Vcut,Vmin-padding,Vmax+padding])
plt.tight_layout()
plt.savefig(plotname)


plt.figure(figsize=(6,5))
if plotvsd==True:
    plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f, Cm=%.2f, l=%i' % (idur,iamp,cm,dendlen))
    plt.plot(time0,V0,label=r'$d=$%s' % str(denddiams[0]))
    for i in range(1,Npl):
        plt.plot(times,Vlist[i],label=r'$d=$%s' % str(denddiams[i]))
elif plotvsl==True:
    plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f, Cm=%.2f, d=%.2f' % (idur,iamp,cm,denddiam))
    plt.plot(time0,V0,label=r'$l=$%s' % str(dendlens[0]))
    for i in range(1,Npl):
        plt.plot(times,Vlist[i],label=r'$l=$%s' % str(dendlens[i]))
plt.xlabel(r'$t$ [ms]')
plt.ylabel(r'$V$ [mV]')
plt.legend(loc='upper right')
plt.axis([idelay-0.2,Vcut_zoom,Vmin-padding,Vmax+padding])
plt.tight_layout()
plt.savefig(plotname_2)

plt.figure(figsize=(6,5))
if plotvsd==True:
    plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f, Cm=%.2f, l=%i' % (idur,iamp,cm,dendlen))
    plt.plot(time0,V0norm,label=r'$d=$%s' % str(denddiams[0]))
    for i in range(1,Npl):
        plt.plot(times,Vnormlist[i],label=r'$d=$%s' % str(denddiams[i]))
    plt.plot(taus,tauheights_norm,'o')
elif plotvsl==True:
    plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f, Cm=%.2f, d=%.2f' % (idur,iamp,cm,denddiam))
    plt.plot(time0,V0norm,label=r'$l=$%s' % str(dendlens[0]))
    for i in range(1,Npl):
        plt.plot(times,Vnormlist[i],label=r'$l=$%s' % str(dendlens[i]))
    plt.plot(taus,tauheights_norm,'o')
plt.xlabel(r'$t$ [ms]')
plt.ylabel(r'$V/V_{max}$')
plt.legend(loc='upper right')
plt.axis([idelay-0.2,Vcut_zoom,0,1.01])
plt.tight_layout()
plt.savefig(plotname_norm)

plt.figure(figsize=(6,5))
if plotvsd==True:
    plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f, Cm=%.2f, l=%i' % (idur,iamp,cm,dendlen))
    plt.plot(time0,V0,label=r'$d=$%s' % str(denddiams[0]))
    plt.plot(timecut,fitlist[0],'--', color='grey')
    for i in range(1,Npl):
        plt.plot(times,Vlist[i],label=r'$d=$%s' % str(denddiams[i]))
        plt.plot(timecut,fitlist[i],'--', color='grey') #label=r'$d=$%s' % str(denddiams[i]))
    plt.plot(taus,tauheights,'o')
elif plotvsl==True:
    plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f, Cm=%.2f, d=%.2f' % (idur,iamp,cm,denddiam))
    plt.plot(time0,V0,label=r'$l=$%s' % str(dendlens[0]))
    plt.plot(timecut,fitlist[0],'--', color='grey')
    for i in range(1,Npl):
        plt.plot(times,Vlist[i],label=r'$l=$%s' % str(dendlens[i]))
        plt.plot(timecut,fitlist[i],'--', color='grey')
    plt.plot(taus,tauheights,'o')
plt.xlabel(r'$t$ [ms]')
plt.ylabel(r'$V$ [mV]')
plt.legend(loc='upper right')
plt.axis([idelay-0.2,Vcut_zoom,Vmin-padding,Vmax+padding])
plt.tight_layout()
plt.savefig(plotname_withfit)

plt.figure(figsize=(6,5)) ## CROSSES
if plotvsd==True:
    plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f, Cm=%.2f, l=%i' % (idur,iamp,cm,dendlen))
    tauplotx = np.zeros(2)
    tauploty = np.zeros(2)
    tauhx    = np.zeros(2)
    tauhy    = np.zeros(2)
    tauhx[0] = idelay-0.2
    tauhx[1] = Vcut_zoom
    tauhy[0] = tauheights[0]
    tauhy[1] = tauheights[0]
    tauplotx[0] = taus[0]
    tauplotx[1] = taus[0]
    tauploty[0] = Vmin-padding
    tauploty[1] = Vmax+padding
    plt.plot(time0,V0,label=r'$d=$%s' % str(denddiams[0]))
    plt.plot(tauplotx,tauploty,'--',color='grey')
    plt.plot(tauhx,tauhy,'--',color='grey')
    for i in range(1,Npl):
        plt.plot(times,Vlist[i],label=r'$d=$%s' % str(denddiams[i]))
        tauhy[0] = tauheights[i]
        tauhy[1] = tauheights[i]
        tauplotx[0] = taus[i]
        tauplotx[1] = taus[i]
        plt.plot(tauplotx,tauploty,'--',color='grey')
        plt.plot(tauhx,tauhy,'--',color='grey')
    plt.plot(taus,tauheights,'o')
elif plotvsl==True:
    plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f, Cm=%.2f, d=%.2f' % (idur,iamp,cm,denddiam))
    tauplotx = np.zeros(2)
    tauploty = np.zeros(2)
    tauhx    = np.zeros(2)
    tauhy    = np.zeros(2)
    tauhx[0] = idelay-0.2
    tauhx[1] = Vcut_zoom
    tauhy[0] = tauheights[0]
    tauhy[1] = tauheights[0]
    tauplotx[0] = taus[0]
    tauplotx[1] = taus[0]
    tauploty[0] = Vmin-padding
    tauploty[1] = Vmax+padding
    plt.plot(time0,V0,label=r'$l=$%s' % str(dendlens[0]))
    plt.plot(tauplotx,tauploty,'--',color='grey')
    plt.plot(tauhx,tauhy,'--',color='grey')
    for i in range(1,Npl):
        plt.plot(times,Vlist[i],label=r'$l=$%s' % str(dendlens[i]))
        tauhy[0] = tauheights[i]
        tauhy[1] = tauheights[i]
        tauplotx[0] = taus[i]
        tauplotx[1] = taus[i]
        plt.plot(tauplotx,tauploty,'--',color='grey')
        plt.plot(tauhx,tauhy,'--',color='grey')
    plt.plot(taus,tauheights,'o')
plt.xlabel(r'$t$ [ms]')
plt.ylabel(r'$V$ [mV]')
plt.legend(loc='upper right')
plt.axis([idelay-0.2,Vcut_zoom,Vmin-padding,Vmax+padding])
plt.tight_layout()
plt.savefig(plotname_crosses)


plt.figure(figsize=(6,5)) ## Shift, only for testing, won't save.
if plotvsd==True:
    plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f, Cm=%.2f, l=%i' % (idur,iamp,cm,dendlen))
    plt.plot(time0,V0shift,label=r'$d=$%s' % str(denddiams[0]))
    for i in range(1,Npl):
        plt.plot(times,Vshiftlist[i],label=r'$d=$%s' % str(denddiams[i]))
elif plotvsl==True:
    plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f, Cm=%.2f, d=%.2f' % (idur,iamp,cm,denddiam))
    plt.plot(time0,V0shift,label=r'$l=$%s' % str(dendlens[0]))
    for i in range(1,Npl):
        plt.plot(times,Vshiftlist[i],label=r'$l=$%s' % str(dendlens[i]))
plt.xlabel(r'$t$ [ms]')
plt.ylabel(r'$V-V_{min}$ [mV]')
plt.legend(loc='upper right')
plt.axis([idelay-0.2,Vcut_zoom,0,Vmax-Vmin+0.5])
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(taus,tauheights,'o')
plt.xlabel(r'$\tau$ (ms)')
plt.ylabel(r'$V$ (mV)')
plt.show()

print('taus:',taus)
print('V0:',V0)