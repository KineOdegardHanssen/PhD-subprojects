import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

# change the default font family
plt.rcParams.update({'font.family':'Arial'})

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

cm         = 1.0
skiptime   = 500
spikedurat = -40
idur       = 1000 #100 # ms
idelay     = 10
v_init     = -86.8 # mV
Ra         = 100
somasize   = 10 # 15 # 
dtexp      = -7

gcahvas = [1.0]
gbks    = [1.0]

for gcahva in gcahvas:
    for gbk in gbks:
        fig = plt.figure(figsize=(18,12),dpi=300)
        
        gs = gridspec.GridSpec(3, 10)
        
        ax1  = plt.subplot(gs[0, 0:2])
        ax2  = plt.subplot(gs[0, 2:4])
        ax3  = plt.subplot(gs[0, 4:6])
        ax4  = plt.subplot(gs[0, 6:8])
        ax5  = plt.subplot(gs[0, 8:10])
        ax6  = plt.subplot(gs[1, 0:2])
        ax7  = plt.subplot(gs[1, 2:4])
        ax8  = plt.subplot(gs[1, 4:6])
        ax9  = plt.subplot(gs[1, 6:8])
        ax10 = plt.subplot(gs[1, 8:10])
        ax11 = plt.subplot(gs[2, 0:2])
        ax12 = plt.subplot(gs[2, 2:4])
        ax13 = plt.subplot(gs[2, 4:6])
        ax14 = plt.subplot(gs[2, 6:8])
        ax15 = plt.subplot(gs[2, 8:10])
        
        #fig.suptitle(r'Properties',fontsize=20)
        
        ax1.set_title(r'A',loc='left',pad=50,fontsize=18)
        ax2.set_title(r'B',loc='left',pad=50,fontsize=18)
        ax3.set_title(r'C',loc='left',pad=50,fontsize=18)
        ax4.set_title(r'D',loc='left',pad=50,fontsize=18)
        ax5.set_title(r'E',loc='left',pad=50,fontsize=18)
        ax6.set_title(r'F',loc='left',pad=50,fontsize=18)
        ax7.set_title(r'G',loc='left',pad=50,fontsize=18)
        ax8.set_title(r'H',loc='left',pad=50,fontsize=18)
        ax9.set_title(r'I',loc='left',pad=50,fontsize=18)
        ax10.set_title(r'J',loc='left',pad=50,fontsize=18)
        ax11.set_title(r'K',loc='left',pad=50,fontsize=18)
        ax12.set_title(r'L',loc='left',pad=50,fontsize=18)
        ax13.set_title(r'M',loc='left',pad=50,fontsize=18)
        ax14.set_title(r'N',loc='left',pad=50,fontsize=18)
        ax15.set_title(r'O',loc='left',pad=50,fontsize=18)
        
        outfolder = 'Compare/Soma%i/' % somasize
        plotname = outfolder+'compare_allmodels_fI_all_mainmodelgs_gcahva'+str(gcahva)+'p_gsk'+str(gbk)+'p.png'
        
        #### Subplot 1 ####    
        model_folders = ['','CaHVA_Allen/gCa2.986e-5/','CaHVA_Allen/bk/gCa2.986e-5_gSK0.0028/']
        labels = ['base','CaHVA','CaHVA, bk']
        infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+ '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt']
        Nmodels = len(model_folders)
        
        ax1.set_title(r'CaHVA (Allen); bk (Hjorth)',fontsize=12)
        for i in range(Nmodels):
            Nspikes      = []
            iamp_Nspikes = []
            
            # Get names
            infolder = model_folders[i]+'Results/Soma%i/' % somasize
            infilename_Nspikes = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'#infilenames[i]
            
            # Open files
            infile_Nspikes = open(infilename_Nspikes,'r')
            
            lines_Nspikes = infile_Nspikes.readlines()
            
            for line in lines_Nspikes:
                words = line.split()
                if len(words)>0:
                    iamp_Nspikes.append(float(words[0]))
                    Nspikes.append(float(words[1]))
            
            ax1.plot(iamp_Nspikes,Nspikes,label=labels[i])
        
            infile_Nspikes.close()
    
        #### Subplot 2 ####
        model_folders = ['','CaHVA_Allen/gCa2.986e-5/','CaHVA_Allen/BK_I_Zhang/gCa2.986e-5_gSK0.0028/']
        labels = ['base','CaHVA','CaHVA, BK']
        infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt']
        Nmodels = len(model_folders)
        
        cm = 1.0
        
        ax2.set_title(r'CaHVA (Allen); BK (Zhang)',fontsize=12)
        for i in range(Nmodels):
            Nspikes      = []
            iamp_Nspikes = []
            
            # Get names
            infolder = model_folders[i]+'Results/Soma%i/' % somasize
            infilename_Nspikes = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'#infilenames[i]
            
            # Open files
            infile_Nspikes = open(infilename_Nspikes,'r')
            
            lines_Nspikes = infile_Nspikes.readlines()
            
            for line in lines_Nspikes:
                words = line.split()
                if len(words)>0:
                    iamp_Nspikes.append(float(words[0]))
                    Nspikes.append(float(words[1]))
            
            ax2.plot(iamp_Nspikes,Nspikes,label=labels[i])
            
            infile_Nspikes.close()
        
        #### Subplot 3 ####
        model_folders = ['','CaHVA_Allen/gCa2.986e-5/','CaHVA_Allen/BK_AitOuares/gCa2.986e-5_gSK0.0028/']
        labels = ['base','CaHVA','CaHVA, BK']
        infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt']
        Nmodels = len(model_folders)
        
        cm = 1.0
        
        ax3.set_title(r'CaHVA (Allen); BK (Ait Ouares)',fontsize=12)
        for i in range(Nmodels):
            Nspikes      = []
            iamp_Nspikes = []
            
            # Get names
            infolder = model_folders[i]+'Results/Soma%i/' % somasize
            infilename_Nspikes = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'#infilenames[i]
            
            # Open files
            infile_Nspikes = open(infilename_Nspikes,'r')
            
            lines_Nspikes = infile_Nspikes.readlines()
            
            for line in lines_Nspikes:
                words = line.split()
                if len(words)>0:
                    iamp_Nspikes.append(float(words[0]))
                    Nspikes.append(float(words[1]))
            
            ax3.plot(iamp_Nspikes,Nspikes,label=labels[i])
            
            infile_Nspikes.close()
        
        
        #### Subplot 4 ####
        model_folders = ['','CaHVA_Allen/gCa2.986e-5/','CaHVA_Allen/SK_AitOuares/gCa2.986e-5_gSK0.0028/']
        labels = ['base','CaHVA','CaHVA, SK']
        infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm) + '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gsk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt']
        Nmodels = len(model_folders)
        
        cm = 1.0
        
        ax4.set_title(r'CaHVA (Allen); SK (Ait Ouares)',fontsize=12)
        for i in range(Nmodels):
            Nspikes      = []
            iamp_Nspikes = []
            
            # Get names
            infolder = model_folders[i]+'Results/Soma%i/' % somasize
            infilename_Nspikes = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'#infilenames[i]
            
            # Open files
            infile_Nspikes = open(infilename_Nspikes,'r')
            
            lines_Nspikes = infile_Nspikes.readlines()
            
            for line in lines_Nspikes:
                words = line.split()
                if len(words)>0:
                    iamp_Nspikes.append(float(words[0]))
                    Nspikes.append(float(words[1]))
            
            ax4.plot(iamp_Nspikes,Nspikes,label=labels[i])
            
            infile_Nspikes.close()

        #### Subplot 5 ####
        model_folders = ['','CaHVA_Allen/gCa2.986e-5/','CaHVA_Allen/SK_Allen/gCa2.986e-5_gSK0.0028/']
        labels = ['base','CaHVA','CaHVA, SK']
        infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gsk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt']
        Nmodels = len(model_folders)
        
        cm = 1.0
        
        ax5.set_title(r'CaHVA (Allen); SK (Allen)',fontsize=12)
        for i in range(Nmodels):
            Nspikes      = []
            iamp_Nspikes = []
            
            # Get names
            infolder = model_folders[i]+'Results/Soma%i/' % somasize
            infilename_Nspikes = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'#infilenames[i]
            
            # Open files
            infile_Nspikes = open(infilename_Nspikes,'r')
            
            lines_Nspikes = infile_Nspikes.readlines()
            
            for line in lines_Nspikes:
                words = line.split()
                if len(words)>0:
                    iamp_Nspikes.append(float(words[0]))
                    Nspikes.append(float(words[1]))
            
            ax5.plot(iamp_Nspikes,Nspikes,label=labels[i])
            
            infile_Nspikes.close()
            
            
        #### Subplot 6 ####
        model_folders = ['','canin_Konstantoudaki/gCa2.986e-5/','canin_Konstantoudaki/bk/gCa2.986e-5_gSK0.0028/']
        labels = ['base','canin','canin, bk']
        infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt']
        Nmodels = len(model_folders)
        
        cm = 1.0
        
        ax6.set_title(r'canin (Konstantoudaki); bk (Hjorth)',fontsize=12)
        for i in range(Nmodels):
            Nspikes      = []
            iamp_Nspikes = []
            
            # Get names
            infolder = model_folders[i]+'Results/Soma%i/' % somasize
            infilename_Nspikes = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'#infilenames[i]
            
            # Open files
            infile_Nspikes = open(infilename_Nspikes,'r')
            
            lines_Nspikes = infile_Nspikes.readlines()
            
            for line in lines_Nspikes:
                words = line.split()
                if len(words)>0:
                    iamp_Nspikes.append(float(words[0]))
                    Nspikes.append(float(words[1]))
            
            ax6.plot(iamp_Nspikes,Nspikes,label=labels[i])
            
            infile_Nspikes.close()
            
            
        #### Subplot 7 ####
        model_folders = ['','canin_Konstantoudaki/gCa2.986e-5/','canin_Konstantoudaki/BK_I_Zhang/gCa2.986e-5_gSK0.0028/']
        labels = ['base','canin','canin, BK']
        infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt']
        Nmodels = len(model_folders)
        
        cm = 1.0
        
        ax7.set_title(r'canin (Konstantoudaki); BK (Zhang)',fontsize=12)
        for i in range(Nmodels):
            Nspikes      = []
            iamp_Nspikes = []
            
            # Get names
            infolder = model_folders[i]+'Results/Soma%i/' % somasize
            infilename_Nspikes = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'#infilenames[i]
            
            # Open files
            infile_Nspikes = open(infilename_Nspikes,'r')
            
            lines_Nspikes = infile_Nspikes.readlines()
            
            for line in lines_Nspikes:
                words = line.split()
                if len(words)>0:
                    iamp_Nspikes.append(float(words[0]))
                    Nspikes.append(float(words[1]))
            
            ax7.plot(iamp_Nspikes,Nspikes,label=labels[i])
            
            infile_Nspikes.close()
            
        
        #### Subplot 8 ####
        model_folders = ['','canin_Konstantoudaki/gCa2.986e-5/','canin_Konstantoudaki/BK_AitOuares/gCa2.986e-5_gSK0.0028/']
        labels = ['base','canin','canin, BK']
        infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt']
        Nmodels = len(model_folders)
        
        cm = 1.0
        
        ax8.set_title(r'canin (Konstantoudaki); BK (Ait Ouares)',loc='right',fontsize=12)
        for i in range(Nmodels):
            Nspikes      = []
            iamp_Nspikes = []
            
            # Get names
            infolder = model_folders[i]+'Results/Soma%i/' % somasize
            infilename_Nspikes = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'#infilenames[i]
            
            # Open files
            infile_Nspikes = open(infilename_Nspikes,'r')
            
            lines_Nspikes = infile_Nspikes.readlines()
            
            for line in lines_Nspikes:
                words = line.split()
                if len(words)>0:
                    iamp_Nspikes.append(float(words[0]))
                    Nspikes.append(float(words[1]))
            
            ax8.plot(iamp_Nspikes,Nspikes,label=labels[i])
            
            infile_Nspikes.close()
        
        #### Subplot 9 ####
        model_folders = ['','canin_Konstantoudaki/gCa2.986e-5/','canin_Konstantoudaki/SK_AitOuares/gCa2.986e-5_gSK0.0028/']
        labels = ['base','canin','canin, BK']
        infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gsk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt']
        Nmodels = len(model_folders)
        
        cm = 1.0
        
        ax9.set_title(r'canin (Konstantoudaki); SK (Ait Ouares)',loc='right',fontsize=12) #'_gsk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p_
        for i in range(Nmodels):
            Nspikes      = []
            iamp_Nspikes = []
            
            # Get names
            infolder = model_folders[i]+'Results/Soma%i/' % somasize
            infilename_Nspikes = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'#infilenames[i]
            
            # Open files
            infile_Nspikes = open(infilename_Nspikes,'r')
            
            lines_Nspikes = infile_Nspikes.readlines()
            
            for line in lines_Nspikes:
                words = line.split()
                if len(words)>0:
                    iamp_Nspikes.append(float(words[0]))
                    Nspikes.append(float(words[1]))
            
            ax9.plot(iamp_Nspikes,Nspikes,label=labels[i])
            
            infile_Nspikes.close()
        
            
        #### Subplot 10 ####
        model_folders = ['','canin_Konstantoudaki/gCa2.986e-5/','canin_Konstantoudaki/SK_Allen/gCa2.986e-5_gSK0.0028/']
        labels = ['base','canin','canin, SK']
        infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gsk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt']
        Nmodels = len(model_folders)
        
        cm = 1.0
        
        ax10.set_title(r'canin (Konstantoudaki); SK (Allen)',loc='right',fontsize=12)
        for i in range(Nmodels):
            Nspikes      = []
            iamp_Nspikes = []
            
            # Get names
            infolder = model_folders[i]+'Results/Soma%i/' % somasize
            infilename_Nspikes = infolder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'#infilenames[i]
            
            # Open files
            infile_Nspikes = open(infilename_Nspikes,'r')
            
            lines_Nspikes = infile_Nspikes.readlines()
            
            for line in lines_Nspikes:
                words = line.split()
                if len(words)>0:
                    iamp_Nspikes.append(float(words[0]))
                    Nspikes.append(float(words[1]))
            
            ax10.plot(iamp_Nspikes,Nspikes,label=labels[i])
            
            infile_Nspikes.close()   
        #### Subplot 11 ####
        model_folders = ['','caq/gCa2.986e-5/','caq/bk/gCa2.986e-5_gSK0.0028/']
        labels = ['base','caq','caq, bk']
        infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt']
        Nmodels = len(model_folders)
        
        cm = 1.0
        
        ax11.set_title(r'caq (Hjorth); bk (Hjorth)',fontsize=10)
        for i in range(Nmodels):
            Nspikes      = []
            iamp_Nspikes = []
            
            # Get names
            infolder = model_folders[i]+'Results/Soma%i/' % somasize
            infilename_Nspikes = infolder+infilenames[i]
            
            # Open files
            infile_Nspikes = open(infilename_Nspikes,'r')
            
            lines_Nspikes = infile_Nspikes.readlines()
            
            for line in lines_Nspikes:
                words = line.split()
                if len(words)>0:
                    iamp_Nspikes.append(float(words[0]))
                    Nspikes.append(float(words[1]))
            
            ax11.plot(iamp_Nspikes,Nspikes,label=labels[i])
            
            infile_Nspikes.close()
            
            
        #### Subplot 12 ####
        model_folders = ['','caq/gCa2.986e-5/','caq/BK_I_Zhang/gCa2.986e-5_gSK0.0028/']
        labels = ['base','caq','caq, BK']
        infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk1.0p' + '_gCaHVA1.0p__Nspikes_vs_I_s'+str(skiptime)+'.txt']
        Nmodels = len(model_folders)
        
        cm = 1.0
        
        ax12.set_title(r'caq (Hjorth); BK (Zhang)',fontsize=10)
        for i in range(Nmodels):
            Nspikes      = []
            iamp_Nspikes = []
            
            # Get names
            infolder = model_folders[i]+'Results/Soma%i/' % somasize
            infilename_Nspikes = infolder+infilenames[i]
            
            # Open files
            infile_Nspikes = open(infilename_Nspikes,'r')
            
            lines_Nspikes = infile_Nspikes.readlines()
            
            for line in lines_Nspikes:
                words = line.split()
                if len(words)>0:
                    iamp_Nspikes.append(float(words[0]))
                    Nspikes.append(float(words[1]))
            
            ax12.plot(iamp_Nspikes,Nspikes,label=labels[i])
            
            infile_Nspikes.close()
            
        
        #### Subplot 13 ####
        model_folders = ['','caq/gCa2.986e-5/','caq/BK_AitOuares/gCa2.986e-5_gSK0.0028/']
        labels = ['base','caq','caq, BK']
        infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt']
        Nmodels = len(model_folders)
        
        cm = 1.0
        
        ax13.set_title(r'caq (Hjorth); BK (Ait Ouares)',loc='right',fontsize=10)
        for i in range(Nmodels):
            Nspikes      = []
            iamp_Nspikes = []
            
            # Get names
            infolder = model_folders[i]+'Results/Soma%i/' % somasize
            infilename_Nspikes = infolder+infilenames[i]
            
            # Open files
            infile_Nspikes = open(infilename_Nspikes,'r')
            
            lines_Nspikes = infile_Nspikes.readlines()
            
            for line in lines_Nspikes:
                words = line.split()
                if len(words)>0:
                    iamp_Nspikes.append(float(words[0]))
                    Nspikes.append(float(words[1]))
            
            ax13.plot(iamp_Nspikes,Nspikes,label=labels[i])
            
            infile_Nspikes.close()
        
        #### Subplot 14 ####
        model_folders = ['','caq/gCa2.986e-5/','caq/SK_AitOuares/gCa2.986e-5_gSK0.0028/']
        labels = ['base','caq','caq, BK']
        infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt']
        Nmodels = len(model_folders)
        
        cm = 1.0
        
        ax14.set_title(r'caq (Hjorth); SK (Ait Ouares)',loc='right',fontsize=10) #'_gsk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p_
        for i in range(Nmodels):
            Nspikes      = []
            iamp_Nspikes = []
            
            # Get names
            infolder = model_folders[i]+'Results/Soma%i/' % somasize
            infilename_Nspikes = infolder+infilenames[i]
            
            # Open files
            infile_Nspikes = open(infilename_Nspikes,'r')
            
            lines_Nspikes = infile_Nspikes.readlines()
            
            for line in lines_Nspikes:
                words = line.split()
                if len(words)>0:
                    iamp_Nspikes.append(float(words[0]))
                    Nspikes.append(float(words[1]))
            
            ax14.plot(iamp_Nspikes,Nspikes,label=labels[i])
            
            infile_Nspikes.close()
        
            
        #### Subplot 15 ####
        model_folders = ['','caq/gCa2.986e-5/','caq/SK_Allen/gCa2.986e-5_gSK0.0028/']
        labels = ['base','caq','caq, SK']
        infilenames = ['somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt','somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_gbk'+str(gbk)+'p' + '_gCaHVA'+str(gcahva)+'p__Nspikes_vs_I_s'+str(skiptime)+'.txt']
        Nmodels = len(model_folders)
        
        cm = 1.0
        
        ax15.set_title(r'caq (Hjorth); SK (Allen)',fontsize=10)
        for i in range(Nmodels):
            Nspikes      = []
            iamp_Nspikes = []
            
            # Get names
            infolder = model_folders[i]+'Results/Soma%i/' % somasize
            infilename_Nspikes = infolder+infilenames[i]
            
            # Open files
            infile_Nspikes = open(infilename_Nspikes,'r')
            
            lines_Nspikes = infile_Nspikes.readlines()
            
            for line in lines_Nspikes:
                words = line.split()
                if len(words)>0:
                    iamp_Nspikes.append(float(words[0]))
                    Nspikes.append(float(words[1]))
            
            ax15.plot(iamp_Nspikes,Nspikes,label=labels[i])
            
            infile_Nspikes.close()
        
        
        ax1.set_xlabel('$I$ (nA)',fontsize=14)
        ax1.set_ylabel('$f$ (Hz)',fontsize=14)
        ax1.legend(loc='upper left')#,ncol=1)
        
        ax2.set_xlabel('$I$ (nA)',fontsize=14)
        ax2.set_ylabel('$f$ (Hz)',fontsize=14)
        ax2.legend(loc='upper left')#,ncol=1)
        
        ax3.set_xlabel('$I$ (nA)',fontsize=14)
        ax3.set_ylabel('$f$ (Hz)',fontsize=14)
        ax3.legend(loc='upper left')#,ncol=1)
        
        ax4.set_xlabel('$I$ (nA)',fontsize=14)
        ax4.set_ylabel('$f$ (Hz)',fontsize=14)
        ax4.legend(loc='upper left')#,ncol=1)
        
        ax5.set_xlabel('$I$ (nA)',fontsize=14)
        ax5.set_ylabel('$f$ (Hz)',fontsize=14)
        ax5.legend(loc='upper left')#,ncol=1)
        
        ax6.set_xlabel('$I$ (nA)',fontsize=14)
        ax6.set_ylabel('$f$ (Hz)',fontsize=14)
        ax6.legend(loc='upper left')#,ncol=1)
        
        ax7.set_xlabel('$I$ (nA)',fontsize=14)
        ax7.set_ylabel('$f$ (Hz)',fontsize=14)
        ax7.legend(loc='upper left')#,ncol=1)
        
        ax8.set_xlabel('$I$ (nA)',fontsize=14)
        ax8.set_ylabel('$f$ (Hz)',fontsize=14)
        ax8.legend(loc='upper left')#,ncol=1)
        
        ax9.set_xlabel('$I$ (nA)',fontsize=14)
        ax9.set_ylabel('$f$ (Hz)',fontsize=14)
        ax9.legend(loc='upper left')#,ncol=1)
        
        ax10.set_xlabel('$I$ (nA)',fontsize=14)
        ax10.set_ylabel('$f$ (Hz)',fontsize=14)
        ax10.legend(loc='upper left')#,ncol=1)
        
        ax11.set_xlabel('$I$ (nA)',fontsize=14)
        ax11.set_ylabel('$f$ (Hz)',fontsize=14)
        ax11.legend(loc='upper left')#,ncol=1)
        
        ax12.set_xlabel('$I$ (nA)',fontsize=14)
        ax12.set_ylabel('$f$ (Hz)',fontsize=14)
        ax12.legend(loc='upper left')#,ncol=1)
        
        ax13.set_xlabel('$I$ (nA)',fontsize=14)
        ax13.set_ylabel('$f$ (Hz)',fontsize=14)
        ax13.legend(loc='upper left')#,ncol=1)
        
        ax14.set_xlabel('$I$ (nA)',fontsize=14)
        ax14.set_ylabel('$f$ (Hz)',fontsize=14)
        ax14.legend(loc='upper left')#,ncol=1)
        
        ax15.set_xlabel('$I$ (nA)',fontsize=14)
        ax15.set_ylabel('$f$ (Hz)',fontsize=14)
        ax15.legend(loc='upper left')#,ncol=1)
        
        
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        plt.savefig(plotname)
#plt.show()
