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

## Flags ##
use_flux	= True			# Running on flux or sophie?
get_los		= True			# Upload Line of Sight Data as well?
write_loc	= 'ss_m1_run1'		# Which directory within /stack_data/ to upload data from

if use_flux == True:
	root = str('/nfs/christoq_ls')
else:
	root = str('/n/Christoq1')


## Functions ##

def ss_recover():
	"""
	This function uploads the pickle files from directory stack_data and configures them into multi dimensional arrays.
	It is meant to work with the self-stacked ensembles.
	"""

	# For now, the "varib" dictionary and "HALODATA" array only need to be uploaded from one halo, say the first.
	# However, we must loop over the ensemble data to combine the stack_data arrays.

	## Load Data ##
	# Make wanted variables global
	global varib,HaloID,Halo_P,Halo_V,Gal_P,Gal_V,Gal_Mags,M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,x_range
	global ENS_CAUMASS,ENS_CAUSURF,ENS_NFWSURF,LOS_CAUMASS,LOS_CAUSURF,LOS_NFWSURF
	global ENS_R,ENS_V,ENS_M,LOS_R,LOS_V,LOS_M,ENS_HVD,LOS_HVD

	global ens_r,ens_caumass,los_r,los_caumass

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
	
	global ens_mbias,ens_mscat,los_mbias,los_mscat,ens_vbias,ens_vscat,los_vbias,los_vscat
	global maLOS_CAUMASS,maLOS_HVD,ens_mfrac,ens_hvdfrac,los_mfrac,los_hvdfrac

	# Defined a Masked array for sometimes zero terms
	maLOS_CAUMASS = ma.masked_array(LOS_CAUMASS,mask=LOS_CAUMASS==0)	# Mask '0' values
	maLOS_HVD = ma.masked_array(LOS_HVD,mask=LOS_HVD==0)			# Mask '0' values

	# Ensemble Mass Fraction Arrays after taking logarithm. 
	ens_mfrac = ma.log(ENS_CAUMASS/M_crit200)
	ens_hvdfrac = ma.log(ENS_HVD/HVD)

	# LOS Mass Fraction Arrays
	array_size = halo_num		# halo_num for horizontal avg first, line_num for vertical avg first. See sites page
	los_mfrac,los_hvdfrac = np.zeros(array_size),np.zeros(array_size)
	for a in range(array_size):
		try:
			los_mfrac[a] = astStats.biweightLocation( ma.copy( ma.log( maLOS_CAUMASS[a]/M_crit200[a]) ), 6.0 )
			los_hvdfrac[a] = astStats.biweightLocation( ma.copy( ma.log( maLOS_HVD[a]/HVD[a]) ), 6.0 )
		except:
			los_mfrac[a] = ma.mean( ma.log( maLOS_CAUMASS[a]/M_crit200[a]) )
			los_hvdfrac[a] = ma.mean( ma.log( maLOS_HVD[a]/HVD[a]) )

	# Bias and Scatter Calculation
	ens_mbias,ens_mscat	= astStats.biweightLocation(ma.copy(ens_mfrac),6.0),astStats.biweightScale(ma.copy(ens_mfrac),9.0)
	los_mbias,los_mscat	= astStats.biweightLocation(ma.copy(los_mfrac),6.0),astStats.biweightScale(ma.copy(los_mfrac),9.0)
	ens_vbias,ens_vscat	= astStats.biweightLocation(ma.copy(ens_hvdfrac),6.0),astStats.biweightScale(ma.copy(ens_hvdfrac),9.0)
	los_vbias,los_vscat	= astStats.biweightLocation(ma.copy(los_hvdfrac),6.0),astStats.biweightScale(ma.copy(los_hvdfrac),9.0)


	return













