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
	global varib,HaloID,Halo_P,Halo_V,Gal_P,Gal_V,Gal_Mags,HaloData,x_range
	global ENS_CAUMASS,ENS_CAUSURF,ENS_NFWSURF,LOS_CAUMASS,LOS_CAUSURF,LOS_NFWSURF
	global ENS_R,ENS_V,ENS_M,LOS_R,LOS_V,LOS_M,ENS_HVD,LOS_HVD

	# Create them as lists
	ENS_CAUMASS,ENS_CAUSURF,ENS_NFWSURF,LOS_CAUMASS,LOS_CAUSURF,LOS_NFWSURF = [],[],[],[],[],[]
	ENS_R,ENS_V,ENS_M,LOS_R,LOS_V,LOS_M,ENS_HVD,LOS_HVD = [],[],[],[],[],[],[],[]

	# Initialization step 
	j = 0
	pkl_file = open(root+'/nkern/Stacking/stack_data/'+write_loc+'/Ensemble_'+str(j)+'_Data.pkl','rb')
	input = pkl.Unpickler(pkl_file)

	varib 							= input.load()
	HaloID,Halo_P,Halo_V,Gal_P,Gal_V,Gal_Mags,HaloData	= input.load()
	stack_data						= input.load()

	ens_r,ens_v,ens_m,ens_hvd,ens_caumass,ens_causurf,ens_nfwsurf,los_r,los_v,los_m,los_hvd,los_caumass,los_causurf,los_nfwsurf,x_range = stack_data
	M_crit200,R_crit200,Z,SRAD,ESRAD,HVD = HaloData	

	# Append stack_data to major lists
	ENS_R.append(ens_r)
	ENS_V.append(ens_v)
	ENS_M.append(ens_m)
	ENS_HVD.append(ens_hvd)
	ENS_CAUMASS.append(ens_caumass)
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
	for i in np.arange(1,100):
		j += 1

		pkl_file = open(root+'/nkern/Stacking/stack_data/'+write_loc+'/Ensemble_'+str(j)+'_Data.pkl','rb')
		input = pkl.Unpickler(pkl_file)

		data1 							= input.load()
		data2							= input.load()
		stack_data						= input.load()

		ens_r,ens_v,ens_m,ens_hvd,ens_caumass,ens_causurf,ens_nfwsurf,los_r,los_v,los_m,los_hvd,los_caumass,los_causurf,los_nfwsurf,x_range = stack_data

		# Append stack_data to major lists
		ENS_R.append(ens_r)
		ENS_V.append(ens_v)
		ENS_M.append(ens_m)
		ENS_HVD.append(ens_hvd)
		ENS_CAUMASS.append(ens_caumass)
		ENS_CAUSURF.append(ens_causurf)
		ENS_NFWSURF.append(ens_nfwsurf)
		LOS_R.append(los_r)
		LOS_V.append(los_v)
		LOS_M.append(los_m)
		LOS_HVD.append(los_hvd)
		LOS_CAUMASS.append(los_caumass)
		LOS_CAUSURF.append(los_causurf)
		LOS_NFWSURF.append(los_nfwsurf)

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
	return













