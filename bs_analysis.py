
## To check for completeness of runs
import numpy as np
for (j,k) in zip(np.arange(1,50),[2,5,10,15,25,50,100]*7):
	for i in range(2100/k):
		try:
			f = open('bs_m0_run'+str(j)+'/Ensemble_'+str(i)+'_Data.pkl','rb')
		except:
			print 'j =',j,',i =',i


# Make M200 of bin using median
line_num = varib['line_num']
gal_num = varib['gal_num']
cell_num = varib['cell_num']
BIN_M200,BIN_HVD = [],[]
for i in range(2100/line_num):
	BIN_M200.append(np.median(M_crit200[i*line_num:(i+1)*line_num]))
	BIN_HVD.append(np.median(HVD[i*line_num:(i+1)*line_num]))

BIN_M200 = np.array(BIN_M200)
BIN_HVD = np.array(BIN_HVD)

# Calc Statistics
ens_mfrac = np.log(ENS_CAUMASS/BIN_M200)
ens_vfrac = np.log(ENS_HVD/BIN_HVD)
ens_mbias = np.around( astStats.biweightLocation(ens_mfrac,6.0), 3) * 100
ens_mscat = np.around( astStats.biweightScale(ens_mfrac,9.0), 3) * 100
ens_vbias = np.around( astStats.biweightLocation(ens_vfrac,6.0), 3) * 100
ens_vscat = np.around( astStats.biweightScale(ens_vfrac,9.0),3) * 100

# Make Mass 1-1 plot
mp.loglog([BIN_M200[0],BIN_M200[-1]],[BIN_M200[0],BIN_M200[-1]],'b')
mp.loglog(BIN_M200,ENS_CAUMASS,'ko',alpha=.7)
mp.xlim(6e13, 1.3e15)
mp.ylim(3.5e13,1.3e15)
mp.xlabel('Median of Bin, Binned on Table Mass',fontsize=15)
mp.ylabel('Ensemble Mass Estimate',fontsize=15)
mp.title('1-1 Correlation, Ngal='+str(gal_num)+', LOS='+str(line_num)+', ENS Richness='+str(line_num*gal_num))
mp.figtext(.7,.2,'Cell #'+str(cell_num)+'\nMethod #'+str(varib['method_num'])+'\nBias = '+str(ens_mbias)+'%\nScatter = '+str(ens_mscat)+'%')
mp.show()

# Make HVD 1-1 plot
mp.loglog([BIN_HVD[0],BIN_HVD[-1]],[BIN_HVD[0],BIN_HVD[-1]],'b')
mp.loglog(BIN_HVD,ENS_HVD,'ko')
mp.xlim(BIN_HVD[-1]-50,BIN_HVD[0]+100)
mp.ylim(BIN_HVD[-1]-50,BIN_HVD[0]+100)
mp.xlabel('Median of Bin, Binned on Table HVD',fontsize=15)
mp.ylabel('Ensemble HVD Estimate',fontsize=15)
mp.title('1-1 Correlation, Ngal='+str(gal_num)+', LOS='+str(line_num)+', ENS Richness='+str(line_num*gal_num))
mp.figtext(.7,.2,'Cell #'+str(cell_num)+'\nMethod #'+str(varib['method_num'])+'\nBias = '+str(ens_vbias)+'%\nScatter = '+str(ens_vscat)+'%')
mp.show()









