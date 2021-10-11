import numpy as np
import matplotlib.pyplot as plt
import math

dtexp = -8

plotvsC  = False # True # 
plotvsE  = True # False # 
plotvsg  = False # True # 

varymech = 'Na' # 'K' # 'leak'
varyE = [-10,20,50,60]
varyg = 'None'
    
varylist = [] # Should be redundant
plotstring = ''
if varyE_bool==True and plotvsC==False:
    varylist = varyE
    plotstring = plotstring + 'E'
elif varyE_bool==False and plotvsC==False:
    varylist = varyg
    plotstring = plotstring + 'g'
else:
    plotstring = plotstring + 'C'
Nvary    = len(varylist)

if varymech=='Na':
    folderstring = 'VaryNa/' 
    plotstring = plotstring + ' Na'
elif varymech=='leak':
    folderstring = 'VaryLeak/'
    plotstring = plotstring + ' leak'
elif varymech=='K':
    folderstring = 'VaryK/'
    plotstring = plotstring + ' K'

# Values of membrane capacitance (parameters):
if plotvsC==True:
    cms = [0.5,0.75,1.0,1.25,1.5]
    Npl = len(cms)
    varylist = cms
    # To choose from: [0.01,0.1,0.5,0.8,1.0,1.2,1.5,2.0,3.0,5.0,7.0,10.0,15.0]
elif plotvsE==True: # Include somaonly in these? A lot of work
    varylist = varyE
    Npl      = len(varylist)
elif plotvsg==True:
    varylist = varyg
    Npl      = len(varylist)


paddit   = True
zerobool = False # DO NOT TOUCH!
if paddit==True:
    padding = 2
else:
    padding = 0
# loop de loop? # Yes.
somasize = 10
dendlen  = 100
denddiam = 1
cm     = 1.0
Ra     = 100
v_init = -70
gpas   = 0.0003
vpas   = -65
long   = True
soma_v_init = -65
dendlen  = 1000
denddiam = 2

time0 = []

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
    plotname = folder +'bashh_cms_idur%i_iamp' % idur+str(iamp)+'_Ra%i_ena'%Ra+str(ena)+'_ek' +str(ek)+'_el' +str(el)+'_gnabar'+str(gnabar)+'_gkbar'+str(gkbar)+'_gl'+str(gl)+'_vinit'+str(v_init)+'_Vt.png'
    plotname_2 = folder +'bashh_cms_idur%i_iamp' % idur+str(iamp)+'_Ra%i_ena'%Ra+str(ena)+'_ek' +str(ek)+'_el' +str(el)+'_gnabar'+str(gnabar)+'_gkbar'+str(gkbar)+'_gl'+str(gl)+'_vinit'+str(v_init)+'_Vt_zoomed.png'
    plotname_norm = folder +'bashh_cms_idur%i_iamp' % idur+str(iamp)+'_Ra%i_ena'%Ra+str(ena)+'_ek' +str(ek)+'_el' +str(el)+'_gnabar'+str(gnabar)+'_gkbar'+str(gkbar)+'_gl'+str(gl)+'_vinit'+str(v_init)+'_Vt_normalized.png'
    plotname_crosses = folder +'bashh_cms_idur%i_iamp' % idur+str(iamp)+'_Ra%i_ena'%Ra+str(ena)+'_ek' +str(ek)+'_el' +str(el)+'_gnabar'+str(gnabar)+'_gkbar'+str(gkbar)+'_gl'+str(gl)+'_vinit'+str(v_init)+'_Vt_withtaus.png'
    plotname_withfit = folder +'bashh_cms_idur%i_iamp' % idur+str(iamp)+'_Ra%i_ena'%Ra+str(ena)+'_ek' +str(ek)+'_el' +str(el)+'_gnabar'+str(gnabar)+'_gkbar'+str(gkbar)+'_gl'+str(gl)+'_vinit'+str(v_init)+'_Vt_withfit.png'
    Vcut = 200
elif plotvsd==True:
    folder = 'Comparemodels/BAS_vs_somaonly_passive/dtexp%i/Varyds/' % dtexp
    plotname = folder +'bashh_cm'+str(cm)+'_varyds_dendlen'+str(dendlen)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i_ena'%Ra+str(ena)+'_ek' +str(ek)+'_el' +str(el)+'_gnabar'+str(gnabar)+'_gkbar'+str(gkbar)+'_gl'+str(gl)+'_vinit'+str(v_init)+'_Vt.png'
    plotname_2 = folder +'bashh_cm'+str(cm)+'_varyds_dendlen'+str(dendlen)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i_ena'%Ra+str(ena)+'_ek' +str(ek)+'_el' +str(el)+'_gnabar'+str(gnabar)+'_gkbar'+str(gkbar)+'_gl'+str(gl)+'_vinit'+str(v_init)+'_Vt_zoomed.png'
    plotname_norm = folder +'bashh_cm'+str(cm)+'_varyds_dendlen'+str(dendlen)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i_ena'%Ra+str(ena)+'_ek' +str(ek)+'_el' +str(el)+'_gnabar'+str(gnabar)+'_gkbar'+str(gkbar)+'_gl'+str(gl)+'_vinit'+str(v_init)+'_Vt_normalized.png'
    plotname_crosses = folder +'bashh_cm'+str(cm)+'_varyds_dendlen'+str(dendlen)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i_ena'%Ra+str(ena)+'_ek' +str(ek)+'_el' +str(el)+'_gnabar'+str(gnabar)+'_gkbar'+str(gkbar)+'_gl'+str(gl)+'_vinit'+str(v_init)+'_Vt_withtaus.png'
    plotname_withfit = folder +'bashh_cm'+str(cm)+'_varyds_dendlen'+str(dendlen)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i_ena'%Ra+str(ena)+'_ek' +str(ek)+'_el' +str(el)+'_gnabar'+str(gnabar)+'_gkbar'+str(gkbar)+'_gl'+str(gl)+'_vinit'+str(v_init)+'_Vt_withfit.png'
    Vcut = idelay+25 # Suitable for cm=1.0. # COULD NEED TO UPDATE THIS. Possibly test for the value of cm.
    Vcut_zoom = idelay+5 # Also need to beware of this
elif plotvsl==True:
    folder = 'Comparemodels/BAS_vs_somaonly_passive/dtexp%i/Varyls/' % dtexp
    plotname = folder +'bashh_cm'+str(cm)+'_varyls_denddiam'+str(denddiam)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i_ena'%Ra+str(ena)+'_ek' +str(ek)+'_el' +str(el)+'_gnabar'+str(gnabar)+'_gkbar'+str(gkbar)+'_gl'+str(gl)+'_vinit'+str(v_init)+'_Vt.png'
    plotname_2 = folder +'bashh_cm'+str(cm)+'_varyls_denddiam'+str(denddiam)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i_ena'%Ra+str(ena)+'_ek' +str(ek)+'_el' +str(el)+'_gnabar'+str(gnabar)+'_gkbar'+str(gkbar)+'_gl'+str(gl)+'_vinit'+str(v_init)+'_Vt_zoomed.png'
    plotname_norm = folder +'bashh_cm'+str(cm)+'_varyls_denddiam'+str(denddiam)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i_ena'%Ra+str(ena)+'_ek' +str(ek)+'_el' +str(el)+'_gnabar'+str(gnabar)+'_gkbar'+str(gkbar)+'_gl'+str(gl)+'_vinit'+str(v_init)+'_Vt_normalized.png'
    plotname_crosses = folder +'bashh_cm'+str(cm)+'_varyls_denddiam'+str(denddiam)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i_ena'%Ra+str(ena)+'_ek' +str(ek)+'_el' +str(el)+'_gnabar'+str(gnabar)+'_gkbar'+str(gkbar)+'_gl'+str(gl)+'_vinit'+str(v_init)+'_Vt_withtaus.png'
    plotname_withfit = folder +'bashh_cm'+str(cm)+'_varyls_denddiam'+str(denddiam)+'_idur%i_iamp' % idur+str(iamp)+'_Ra%i_ena'%Ra+str(ena)+'_ek' +str(ek)+'_el' +str(el)+'_gnabar'+str(gnabar)+'_gkbar'+str(gkbar)+'_gl'+str(gl)+'_vinit'+str(v_init)+'_Vt_withfit.png'
    Vcut = idelay+25 # Suitable for cm=1.0. # COULD NEED TO UPDATE THIS. Possibly test for the value of cm.
    Vcut_zoom = idelay+5 # Also need to beware of this

Vlist   = []
fitlist = []
Vnormlist  = []
Vshiftlist = []
Vnormmin   = 1

for i in range(Npl):
    ### Setting file name according to what you want to vary
    elem = varylist[i]
    changestring =''
    if plotvsC==True:
        cm = cms[i]
        infilename  = folder+'basHH_cmf'+str(cmfac)+'_idur%i_iamp' % idur+str(iamp)+'_Ra'+str(Ra)+'_vinit'+str(v_init)+'_V.txt'
    elif plotvsE==True:
        varyE = elem
        changestring = changestring+'_E'+str(varyE)+'_gdflt'
        folder = 'Ball-and-stick models/BAS_passive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam) +'/current_idur%i_iamp'% idur+str(iamp)+'/'
        infilename  = folder+'basHH_cmf'+str(cmfac)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+changestring+'_V_sprxf.txt' 
    elif plotvsg==True:
        varyg = elem
        changestring = changestring+'_Edefault_g'+str(varyg)
        folder = 'Ball-and-stick models/BAS_passive/Results/IStim/Soma%i/dendlen%i/denddiam'% (somasize,dendlen)+str(denddiam) +'/current_idur%i_iamp'% idur+str(iamp)+'/dtexp%i/' % dtexp
        infilename  = folder+'basHH_cmf'+str(cmfac)+'_idur%i_iamp'%idur+str(iamp)+'_Ra%i_vinit' %Ra+str(v_init)+changestring+'_V_sprxf.txt' 
    
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
elif plotvsE==True:
    plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f, Cm=%.2f, l=%i' % (idur,iamp,cm,dendlen))
    plt.plot(time0,V0,label=r'$d=$%s' % str(varylist[0]))
    for i in range(1,Npl):
        plt.plot(times,Vlist[i],label=r'$E=$%s' % str(varylist[i]))
elif plotvsg==True:
    plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f, Cm=%.2f, d=%.2f' % (idur,iamp,cm,denddiam))
    plt.plot(time0,V0,label=r'$l=$%s' % str(varylist[0]))
    for i in range(1,Npl):
        plt.plot(times,Vlist[i],label=r'$g=$%s' % str(varylist[i]))
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
plt.title(r'$V$ vs $t$, idur=%i, iamp=%.1f, Cm=%.2f, %s' % (idur,iamp,cm,plotstring))
if plotvsd==True:
    plt.plot(time0,V0shift,label=r'$d=$%s' % str(denddiams[0]))
    for i in range(1,Npl):
        plt.plot(times,Vshiftlist[i],label=r'$d=$%s' % str(denddiams[i]))
elif plotvsl==True:
    plt.plot(time0,V0shift,label=r'$l=$%s' % str(dendlens[0]))
    for i in range(1,Npl):
        plt.plot(times,Vshiftlist[i],label=r'$l=$%s' % str(dendlens[i]))
plt.xlabel(r'$t$ [ms]')
plt.ylabel(r'$V$ [mV]')
plt.legend(loc='upper right')
plt.axis([idelay-0.2,Vcut_zoom,0,Vmax-Vmin+0.5])
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(taus,tauheights,'o')
plt.show()

print('taus:',taus)
print('V0:',V0)