import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

# change the default font family
plt.rcParams.update({'font.family':'Arial'})

plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('legend', fontsize=15)

testmodels  = [478513437,488462965,478513407] 
cm         = 1.0
spikedurat = -40
idur       = 2000 #100 # ms
iamps      = [0.7,0.6,0.4] # Change from model to model
idelay     = 100
v_init     = -83.7 # mV
Ra         = 100
somasize   = 10 # 15 # 
dendlen    = 1000
denddiam   = 1
nsegments  = 200 

varyE_bool = False
namestringfirst = ''
varygbool    = True
varyIh       = False # True # 
vary_NaV     = False # True # 
vary_Kd      = False # True # 
vary_Kv2like = False # True # 
vary_Kv3_1   = False # True # 
vary_K_T     = False # True # 
vary_Im_v2   = False # True # False #
vary_SK      = False # True # 
vary_Ca_HVA  = False # True # 
vary_Ca_LVA  = False # True # False
vary_gpas    = True # False # True # 

varyg = [0.1,0.3,0.5,0.7,0.9,1.0,1.1,1.5,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0]
Ng = len(varyg)

outfolder = 'figures/Comparemodels/Varygs/'

def reldiffs_one(fdefault,fother):
    rdout  = 0
    rdout2 = 0
    if fdefault!=0 and fother!=0:
        rdout = ((fdefault-fother)/float(fdefault))
        rdout2 = ((fdefault-fother)/float(fother))
    return rdout, rdout2

varylist = [] # Should be redundant
plotstring = '_vary'
varylist = varyg
plotstring = plotstring + 'g'

if varyIh==True:
    namestringfirst = namestringfirst + '_gIh'
    plotstring = plotstring + '_gIh'
    texstring = r'$\bar{g}_{\mathregular{Ih}}$'
if vary_NaV==True:
    namestringfirst = namestringfirst + '_gNaV'
    plotstring = plotstring + '_gNaV'
    texstring = r'$\bar{g}_{\mathregular{NaV}}$'
if vary_Kd==True:
    namestringfirst = namestringfirst + '_gKd'
    plotstring = plotstring + '_gKd'
    texstring = r'$\bar{g}_{\mathregular{Kd}}$'
if vary_Kv2like==True:
    namestringfirst = namestringfirst + '_gKv2like'
    plotstring = plotstring + '_gKv2like'
    texstring = r'$\bar{g}_{\mathregular{Kv2like}}$'
if vary_Kv3_1==True:
    namestringfirst = namestringfirst + '_gKv31'
    plotstring = plotstring + '_gKv31'
    texstring = r'$\bar{g}_{\mathregular{Kv3.1}}$'
if vary_K_T==True:
    namestringfirst = namestringfirst + '_gKT'
    plotstring = plotstring + '_gKT'
    texstring = r'$\bar{g}_{\mathregular{KT}}$'
if vary_Im_v2==True:
    namestringfirst = namestringfirst + '_gImv2'
    plotstring = plotstring + '_gImv2'
    texstring = r'$\bar{g}_{\mathregular{Im,v2}}$'
if vary_SK==True:
    namestringfirst = namestringfirst + '_gSK'
    plotstring = plotstring + '_gSK'
    texstring = r'$\bar{g}_{\mathregular{SK}}$'
if vary_Ca_HVA==True:
    namestringfirst = namestringfirst + '_gCaHVA'
    plotstring = plotstring + '_gCaHVA'
    texstring = r'$\bar{g}_{\mathregular{CaHVA}}$'
if vary_Ca_LVA==True:
    namestringfirst = namestringfirst + '_gCaLVA'
    plotstring = plotstring + '_gCaLVA'
    texstring = r'$\bar{g}_{\mathregular{CaLVA}}$'
if vary_gpas==True: 
    namestringfirst = namestringfirst + '_gpas'
    plotstring = plotstring + '_gpas'
    texstring = r'$g_{\mathregular{pas}}$'

plotname = outfolder+'changing'+namestringfirst+'.png'

###########################################################
j = 0 
g_all       = []
Nspike_all  = []
defaultinds = []
for testmodel in testmodels:
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
    
    # Will not know number of g's in advance (should have the same no. for each model)
    gs      = []
    Nspikes = []
    
    # Set names
    infolder  = 'figures/%i/current_idur%i_iamp' % (testmodel,idur) + str(iamps[j])+'/'
    infilename = infolder+'%s_%i_cmfall'%(namestringfirst,testmodel)+str(cm)+'_idur%i_varyg'% idur+'_manual_Nspikes_vs_g.txt'
    print('infilename:',infilename)
    infile = open(infilename,'r')
    lines = infile.readlines()
    
    i = 0
    for line in lines:
        words = line.split()
        if len(words)>0:
            thisg = float(words[0])
            print('thisg:',thisg)
            gs.append(thisg)
            Nspikes.append(int(words[1]))
            if thisg==1.0:
                defaultinds.append(i)
        i+=1
    g_all.append(gs)
    Nspike_all.append(Nspikes)
    j+=1

# Finding differences (should double check)

gs_all    = []
diffs_all = []
for k in range(len(testmodels)):
    print('k:',k)
    glist    = g_all[k]
    Nspikes  = Nspike_all[k]
    defind   = defaultinds[k]
    fdefault = Nspikes[defind]
    diffs    = []
    gs       = []
    for i in range(len(glist)):
        g = glist[i]
        if g!=1.0:
            rd1,rd2 = reldiffs_one(fdefault,Nspikes[i])
            diffs.append(rd1)
            gs.append(g)
    gs_all.append(gs)
    diffs_all.append(diffs)

# Unpacking for plotting:
gs1    = gs_all[0]
diffs1 = diffs_all[0]
gs2    = gs_all[1]
diffs2 = diffs_all[1]
gs3    = gs_all[2]
diffs3 = diffs_all[2]

# Making labels for plot:
glist1 = []
glist2 = []
glist3 = []
for g in gs1:
    glist1.append(str(g))
for g in gs2:
    glist2.append(str(g))
for g in gs3:
    glist3.append(str(g))

# Plotting:
barWidth = 0.15

# Model 1
br1 = np.arange(len(glist1))
br2 = [x+0.5*barWidth for x in br1]
brcenter = br1

fig = plt.figure(figsize=(30,15),dpi=300)
gs = gridspec.GridSpec(1, 6)
ax1 = plt.subplot(gs[0, 0:2])
ax1.set_title('A',loc='left',fontsize=18)
ax1.set_title('%s, model %i' % (texstring,testmodels[0]))
ax1.bar(br1,diffs1, width=barWidth)
ax1.axhline(y=0,color='k',linestyle='-',linewidth=1.0)
ax1.set_xlabel('%s' % texstring,fontsize=16)
ax1.set_ylabel('Relative difference at max. current',fontsize=16)
plt.xticks(brcenter, glist1)


# Model 2
br1 = np.arange(len(glist2))
br2 = [x+0.5*barWidth for x in br1]
brcenter = br1

ax2 = plt.subplot(gs[0, 2:4])
ax2.set_title('B',loc='left',fontsize=18)
ax2.set_title('%s, model %i' % (texstring,testmodels[1]))
ax2.bar(br1,diffs2, width=barWidth)
ax2.axhline(y=0,color='k',linestyle='-',linewidth=1.0)
ax2.set_xlabel('%s' % texstring,fontsize=16)
ax2.set_ylabel('Relative difference at max. current',fontsize=16)
plt.xticks(brcenter, glist2)


# Model 3
br1 = np.arange(len(glist3))
br2 = [x+0.5*barWidth for x in br1]
brcenter = br1

ax3 = plt.subplot(gs[0, 4:6])
ax3.set_title('C',loc='left',fontsize=18)
ax3.set_title('%s, model %i' % (texstring,testmodels[2]))
ax3.bar(br1,diffs3, width=barWidth)
ax3.axhline(y=0,color='k',linestyle='-',linewidth=1.0)
ax3.set_xlabel('%s' % texstring,fontsize=16)
ax3.set_ylabel('Relative difference at max. current',fontsize=16)
plt.xticks(brcenter, glist3)

fig.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.savefig(plotname)
plt.show()

    
