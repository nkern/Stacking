"""
This program loads previously written data files from flux batch job running on the stacking method.
The data is serialized using cPickle.
"""

## Import Modules ##
import cPickle as pkl
import numpy as np
import matplotlib.pyplot as mp
import numpy.ma as ma
import astStats
import sys
from AttrDict import AttrDict

## Flags ##
use_flux	= True			# Running on flux or sophie?
get_los		= True			# Upload Line of Sight Data as well?
write_loc	= 'ss_m1_run1'		# Which directory within /stack_data/ to upload data from

if use_flux == True:
	root = str('/nfs/christoq_ls')
else:
	root = str('/n/Christoq1')


## Functions ##
class Recover():
	'''This class contains functions that recover serialized data'''

	def __init__(self):
		pass

	def ss_recover(self,write_loc=write_loc,go_global=True):
		"""
		This function uploads the pickle files from directory stack_data and configures them into multi dimensional arrays.
		It is meant to work with the self-stacked ensembles.
		go_global = True makes variables uploaded to global dictionary, False makes it returned to a dictionary
		write_loc = place where data lives
		"""
		# For now, the "varib" dictionary and "HALODATA" array only need to be uploaded from one halo, say the first.
		# However, we must loop over the ensemble data to combine the stack_data arrays.

		## Load Data ##
		# Make wanted variables global

		# Create them as lists
		ENS_CAUMASS,ENS_CAUSURF,ENS_NFWSURF,LOS_CAUMASS,LOS_CAUSURF,LOS_NFWSURF = [],[],[],[],[],[]
		ENS_R,ENS_V,ENS_M,LOS_R,LOS_V,LOS_M,ENS_HVD,LOS_HVD = [],[],[],[],[],[],[],[]

		# Initialization step 
		pkl_file = open(root+'/nkern/Stacking/stack_data/'+write_loc+'/Ensemble_'+str(0)+'_Data.pkl','rb')
		input = pkl.Unpickler(pkl_file)

		stack_data 						= input.load()
		varib							= input.load()
		HaloID,Halo_P,Halo_V,Gal_P,Gal_V,Gal_Mags,HaloData	= input.load()

		ens_r,ens_v,ens_m,ens_hvd,ens_caumass,ens_causurf,ens_nfwsurf,los_r,los_v,los_m,los_hvd,los_caumass,los_causurf,los_nfwsurf,x_range = stack_data
		M_crit200,R_crit200,Z,SRAD,ESRAD,HVD = HaloData	

		# Append stack_data to major lists
		ENS_R.append(ens_r)
		ENS_V.append(ens_v)
		ENS_M.append(ens_m)
		ENS_HVD.append(float(ens_hvd))
		ENS_CAUMASS.append(float(ens_caumass))
		ENS_CAUSURF.append(ens_causurf)
		ENS_NFWSURF.append(ens_nfwsurf)
		LOS_R.append(los_r)
		LOS_V.append(los_v)
		LOS_M.append(los_m)
		LOS_HVD.append(los_hvd)
		LOS_CAUMASS.append(los_caumass)
		LOS_CAUSURF.append(los_causurf)
		LOS_NFWSURF.append(los_nfwsurf)

		# Loop over ensembles
		j = 0
		for i in np.arange(1,100):
			# Progress Bar
			sys.stdout.write("Progress... "+str(j)+"\r")
			sys.stdout.flush()
			j += 1

			pkl_file = open(root+'/nkern/Stacking/stack_data/'+write_loc+'/Ensemble_'+str(i)+'_Data.pkl','rb')
			input = pkl.Unpickler(pkl_file)

			stack_data = input.load()

			ens_r,ens_v,ens_m,ens_hvd,ens_caumass,ens_causurf,ens_nfwsurf,los_r,los_v,los_m,los_hvd,los_caumass,los_causurf,los_nfwsurf,x_range = stack_data


			# Append stack_data to major lists
			ENS_R.append(ens_r)
			ENS_V.append(ens_v)
			ENS_M.append(ens_m)
			ENS_HVD.append(float(ens_hvd))
			ENS_CAUMASS.append(float(ens_caumass))
			ENS_CAUSURF.append(ens_causurf)
			ENS_NFWSURF.append(ens_nfwsurf)
			LOS_R.append(los_r)
			LOS_V.append(los_v)
			LOS_M.append(los_m)
			LOS_HVD.append(los_hvd)
			LOS_CAUMASS.append(los_caumass)
			LOS_CAUSURF.append(los_causurf)
			LOS_NFWSURF.append(los_nfwsurf)
		
		print ''

		# Convert to arrays
		ENS_R = np.array(ENS_R)
		ENS_V = np.array(ENS_V)
		ENS_M = np.array(ENS_M)
		ENS_HVD = np.array(ENS_HVD)
		ENS_CAUMASS = np.array(ENS_CAUMASS)
		ENS_CAUSURF = np.array(ENS_CAUSURF)
		ENS_NFWSURF = np.array(ENS_NFWSURF)
		LOS_R = np.array(LOS_R)
		LOS_V = np.array(LOS_V)
		LOS_M = np.array(LOS_M)
		LOS_HVD = np.array(LOS_HVD)
		LOS_CAUMASS = np.array(LOS_CAUMASS)
		LOS_CAUSURF = np.array(LOS_CAUSURF)
		LOS_NFWSURF = np.array(LOS_NFWSURF)

		#Update global namespace with varib attributes
		globals().update(varib)

		### Statistical Calculations ###

		# Defined a Masked array for sometimes zero terms
		epsilon = 1.0
		maLOS_CAUMASS = ma.masked_array(LOS_CAUMASS,mask=LOS_CAUMASS<epsilon)	# Mask '0' values
		maLOS_HVD = ma.masked_array(LOS_HVD,mask=LOS_HVD<epsilon)		# Mask '0' values

		# Ensemble Mass Fraction Arrays after taking logarithm. 
		ens_mfrac = ma.log(ENS_CAUMASS/M_crit200)
		ens_vfrac = ma.log(ENS_HVD/HVD)

		# LOS Mass Fraction Arrays
		array_size = halo_num		# halo_num for horizontal avg first, line_num for vertical avg first. See sites page
		los_mfrac,los_vfrac = np.zeros(array_size),np.zeros(array_size)
		for a in range(array_size):
			try:
				los_mfrac[a] = astStats.biweightLocation( ma.copy( ma.log( maLOS_CAUMASS[a]/M_crit200[a]) ), 6.0 )
				los_vfrac[a] = astStats.biweightLocation( ma.copy( ma.log( maLOS_HVD[a]/HVD[a]) ), 6.0 )
			except:
				los_mfrac[a] = ma.mean( ma.log( maLOS_CAUMASS[a]/M_crit200[a]) )
				los_vfrac[a] = ma.mean( ma.log( maLOS_HVD[a]/HVD[a]) )

		# Bias and Scatter Calculation
		ens_mbias,ens_mscat	= astStats.biweightLocation(ma.copy(ens_mfrac),6.0),astStats.biweightScale(ma.copy(ens_mfrac),9.0)
		los_mbias,los_mscat	= astStats.biweightLocation(ma.copy(los_mfrac),6.0),astStats.biweightScale(ma.copy(los_mfrac),9.0)
		ens_vbias,ens_vscat	= astStats.biweightLocation(ma.copy(ens_vfrac),6.0),astStats.biweightScale(ma.copy(ens_vfrac),9.0)
		los_vbias,los_vscat	= astStats.biweightLocation(ma.copy(los_vfrac),6.0),astStats.biweightScale(ma.copy(los_vfrac),9.0)

		# Extra Arrays
		avgLOS_CAUMASS = np.array(map(np.median,maLOS_CAUMASS))

		# Return to Namespaces depending on go_global
		names = ['varib','HaloID','Halo_P','Halo_V','Gal_P','Gal_V','Gal_Mags','M_crit200','R_crit200','Z','SRAD','ESRAD','HVD','x_range','ENS_CAUMASS','ENS_CAUSURF','ENS_NFWSURF','LOS_CAUMASS','LOS_CAUSURF','LOS_NFWSURF','ENS_R','ENS_V','ENS_M','LOS_R','LOS_V','LOS_M','ENS_HVD','LOS_HVD','ens_r','ens_caumass','los_r','los_caumass','ens_mbias','ens_mscat','los_mbias','los_mscat','ens_vbias','ens_vscat','los_vbias','los_vscat','maLOS_CAUMASS','maLOS_HVD','avgLOS_CAUMASS','ens_mfrac','ens_vfrac','los_mfrac','los_vfrac']
		data = [varib,HaloID,Halo_P,Halo_V,Gal_P,Gal_V,Gal_Mags,M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,x_range,ENS_CAUMASS,ENS_CAUSURF,ENS_NFWSURF,LOS_CAUMASS,LOS_CAUSURF,LOS_NFWSURF,ENS_R,ENS_V,ENS_M,LOS_R,LOS_V,LOS_M,ENS_HVD,LOS_HVD,ens_r,ens_caumass,los_r,los_caumass,ens_mbias,ens_mscat,los_mbias,los_mscat,ens_vbias,ens_vscat,los_vbias,los_vscat,maLOS_CAUMASS,maLOS_HVD,avgLOS_CAUMASS,ens_mfrac,ens_vfrac,los_mfrac,los_vfrac]
		mydict = dict( zip(names,data) )
		if go_global == True:
			globals().update(mydict)
			return
		elif go_global == False:
			return  mydict	


	def bs_recover():
		pass



class Work(Recover):
	'''This class contains functions that works with the data previously loaded'''

	def __init__(self,Recover):
		pass	

	def statistic_all(self,iter_array=None):
		'''
		This iterates over different richness geometry configurations and runs statistics on data.
		It is recommended to do any calculations (statistics, plots etc.) within the for loop, and
		then feed results back out via a global variable.
		'''

		# Configure Variables
		if iter_array == None:
			iter_array = np.arange(1,21)

		# Define globals
		global ENS_MBIAS,ENS_MSCAT,ENS_VBIAS,ENS_VSCAT,LOS_MBIAS,LOS_MSCAT,LOS_VBIAS,LOS_VSCAT
		global RUN_NUM,GAL_NUM,LINE_NUM

		# Create arrays
		ENS_MBIAS,ENS_MSCAT,ENS_VBIAS,ENS_VSCAT = [],[],[],[]
		LOS_MBIAS,LOS_MSCAT,LOS_VBIAS,LOS_VSCAT = [],[],[],[]
		RUN_NUM,GAL_NUM,LINE_NUM = [],[],[]

		# Iterate over runs
		print '...Loading Data from '+str(len(iter_array))+' runs'
		for i in iter_array:
			print ''
			print 'Working on Run #'+str(i)
			print '-'*25

			write_loc = 'ss_m1_run'+str(i)
			mydict = self.ss_recover(write_loc=write_loc,go_global=False)
			d = AttrDict(mydict)
			RUN_NUM.append(i)
			GAL_NUM.append(d.varib['gal_num'])
			LINE_NUM.append(d.varib['line_num'])
			ENS_MBIAS.append(d.ens_mbias)
			ENS_MSCAT.append(d.ens_mscat)
			ENS_VBIAS.append(d.ens_vbias)
			ENS_VSCAT.append(d.ens_vscat)
			LOS_MBIAS.append(d.los_mbias)
			LOS_MSCAT.append(d.los_mscat)
			LOS_VBIAS.append(d.los_vbias)
			LOS_VSCAT.append(d.los_vscat)	
			
			del mydict,d


	def oto(self,xarray,yarray,style='o',alpha=None):
		'''Simple log log one to one plot setup'''
		p1, = mp.plot(xarray,yarray,style,alpha=alpha)
		mp.plot([xarray[0],xarray[-1]],[xarray[0],xarray[-1]],'b')
		mp.xscale('log')
		mp.yscale('log')
		return p1


## Initialize Classes
R = Recover()
W = Work(Recover)



