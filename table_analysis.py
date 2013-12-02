import cPickle as pkl
import matplotlib.pyplot as mp
import numpy.ma as ma
import numpy as np

file = open('table_analysis_2124halo.pkl','rb')
input = pkl.Unpickler(file)
data = input.load()

(ENS_MBIAS,ENS_MSCAT,ENS_VBIAS,ENS_VSCAT,LOS_MBIAS,LOS_MSCAT,LOS_VBIAS,LOS_VSCAT,RUN_NUM,GAL_NUM,LINE_NUM,RICH_NUM) = data


## Gifford et al. 2013 Guo 2010 Individual Cluster Stats
IND_RICH 		= np.array([10,15,25,35,50,75,100,150],int)
IND_MSCAT		= np.array([.94,.73,.53,.46,.41,.36,.33,.31])
IND_MBIAS 		= np.array([-.64,-.39,-.16,-.10,-.03,.00,.01,.01])
IND_MBIAS_ERR 		= np.array([.09,.08,.05,.05,.04,.03,.03,.03])
##

raise NameError

mp.errorbar(IND_RICH,IND_MBIAS,yerr=IND_MBIAS_ERR,lw=2,c='DarkBlue',alpha=.7)

## ENSEMBLE MASS SCATTER

fig,axes = mp.subplots(1,1)
p1, = axes.plot(RICH_NUM[0],ENS_MSCAT[0],c='DarkBlue',marker='o',alpha=.8)
axes.plot(RICH_NUM[0],ENS_MSCAT[0],'DarkBlue',alpha=.8)
p2, = axes.plot(RICH_NUM[1],ENS_MSCAT[1],c='DarkGreen',marker='o',alpha=.8)
axes.plot(RICH_NUM[1],ENS_MSCAT[1],'DarkGreen',alpha=.8)
p3, = axes.plot(RICH_NUM[2],ENS_MSCAT[2],c='FireBrick',marker='o',alpha=.8)
axes.plot(RICH_NUM[2],ENS_MSCAT[2],'FireBrick',alpha=.8)
p4, = axes.plot(RICH_NUM[3],ENS_MSCAT[3],c='Indigo',marker='o',alpha=.8)
axes.plot(RICH_NUM[3],ENS_MSCAT[3],'Indigo',alpha=.8)
p5, = axes.plot(RICH_NUM[4],ENS_MSCAT[4],c='DarkSlateGrey',marker='o',alpha=.8)
axes.plot(RICH_NUM[4],ENS_MSCAT[4],'DarkSlateGrey',alpha=.8)
p6, = axes.plot(IND_RICH,IND_MSCAT,'ko',ms=8,alpha=.8)
axes.plot(IND_RICH,IND_MSCAT,'k',lw=3,alpha=.4)
axes.set_xlabel('Ensemble Richness, N, within R200',fontsize=14)
axes.set_ylabel('Mass Estimate Percent Scatter',fontsize=14)
axes.set_title('Self Stacked Ensemble Mass Scatter vs. Richness, 12/2/13')
axes.legend([p1,p2,p3,p4,p5,p6],["N=5, 2<LOS<100","N=10, 2<LOS<100","N=15, 2<LOS<100","N=25, 2<LOS<100","N=50, 2<LOS<100","Gifford 2013"],loc=1)
axes.set_xlim(-10,1000)
axes.set_ylim(0,1.0)
axes.locator_params(axis='both',nbins=12)
##



## ENSEMBLE MASS BIAS
mp.figure()
mp.plot(RICH_NUM[0],ENS_MBIAS[0],c='DarkBlue',marker='o',alpha=.7)
p1, = mp.plot(RICH_NUM[0],ENS_MBIAS[0],'DarkBlue',alpha=.7)
mp.plot(RICH_NUM[1],ENS_MBIAS[1],c='DarkGreen',marker='o',alpha=.7)
p2, = mp.plot(RICH_NUM[1],ENS_MBIAS[1],'DarkGreen',alpha=.7)
mp.plot(RICH_NUM[2],ENS_MBIAS[2],c='FireBrick',marker='o',alpha=.7)
p3, = mp.plot(RICH_NUM[2],ENS_MBIAS[2],'FireBrick',alpha=.7)
mp.plot(RICH_NUM[3],ENS_MBIAS[3],c='Indigo',marker='o',alpha=.7)
p4, = mp.plot(RICH_NUM[3],ENS_MBIAS[3],'Indigo',alpha=.7)
mp.plot(RICH_NUM[4],ENS_MBIAS[4],c='DarkSlateGrey',marker='o',alpha=.7)
p5, = mp.plot(RICH_NUM[4],ENS_MBIAS[4],'DarkSlateGrey',alpha=.7)
mp.xlabel('Ensemble Richness within R200',fontsize=14)
mp.ylabel('Mass Estimate Bias',fontsize=14)
mp.title('Mass Bias vs. Richness, 11/24/13')
mp.legend([p1,p2,p3,p4,p5],["Row1","Row2","Row3","Row4","Row5"],loc=4)
mp.ylim(-.9,.1)
##


## LOS MASS SCATTER
mp.figure()
mp.plot(GAL_NUM[0],LOS_MSCAT[0],c='DarkBlue',marker='o',alpha=.7)
p1, = mp.plot(GAL_NUM[0],LOS_MSCAT[0],'DarkBlue',alpha=.7)
mp.plot(GAL_NUM[1],LOS_MSCAT[1],c='DarkGreen',marker='o',alpha=.7)
p2, = mp.plot(GAL_NUM[1],LOS_MSCAT[1],'DarkGreen',alpha=.7)
mp.plot(GAL_NUM[2],LOS_MSCAT[2],c='FireBrick',marker='o',alpha=.7)
p3, = mp.plot(GAL_NUM[2],LOS_MSCAT[2],'FireBrick',alpha=.7)
mp.plot(GAL_NUM[3],LOS_MSCAT[3],c='Indigo',marker='o',alpha=.7)
p4, = mp.plot(GAL_NUM[3],LOS_MSCAT[3],'Indigo',alpha=.7)
mp.plot(GAL_NUM[4],LOS_MSCAT[4],c='DarkSlateGrey',marker='o',alpha=.7)
p5, = mp.plot(GAL_NUM[4],LOS_MSCAT[4],'DarkSlateGrey',alpha=.7)
mp.xlabel('LOS Total Sample Richness within R200',fontsize=14)
mp.ylabel('Mass Estimate Scatter',fontsize=14)
mp.title('Mass Scatter vs. Richness, 11/24/13')
mp.legend([p1,p2,p3,p4,p5],["Row1","Row2","Row3","Row4","Row5"],loc=1)
mp.ylim(0,1.0)
##


## LOS MASS BIAS
mp.figure()
mp.plot(GAL_NUM[0],LOS_MBIAS[0],c='DarkBlue',marker='o',alpha=.7)
p1, = mp.plot(GAL_NUM[0],LOS_MBIAS[0],'DarkBlue',alpha=.7)
mp.plot(GAL_NUM[1],LOS_MBIAS[1],c='DarkGreen',marker='o',alpha=.7)
p2, = mp.plot(GAL_NUM[1],LOS_MBIAS[1],'DarkGreen',alpha=.7)
mp.plot(GAL_NUM[2],LOS_MBIAS[2],c='FireBrick',marker='o',alpha=.7)
p3, = mp.plot(GAL_NUM[2],LOS_MBIAS[2],'FireBrick',alpha=.7)
mp.plot(GAL_NUM[3],LOS_MBIAS[3],c='Indigo',marker='o',alpha=.7)
p4, = mp.plot(GAL_NUM[3],LOS_MBIAS[3],'Indigo',alpha=.7)
mp.plot(GAL_NUM[4],LOS_MBIAS[4],c='DarkSlateGrey',marker='o',alpha=.7)
p5, = mp.plot(GAL_NUM[4],LOS_MBIAS[4],'DarkSlateGrey',alpha=.7)
mp.xlabel('LOS Total Sample Richness within R200',fontsize=14)
mp.ylabel('Mass Estimate Bias',fontsize=14)
mp.title('Mass Bias vs. Richness, 11/24/13')
mp.legend([p1,p2,p3,p4,p5],["Row1","Row2","Row3","Row4","Row5"],loc=4)
mp.ylim(-.9,.1)
##


mp.show()
















