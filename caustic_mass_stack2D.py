## caustic_mass_stack2D.py
'''
########
-This program uses the Caustic Technique to estimate the mass of galaxy clusters, after having applied a stacking technique to create an ensemble cluster.
-This code works in tandem with caustic.class_stack2D.py
-Parts of the program are either modelled are taken directly from D. Gifford's CausticMass.py code.
-This is the most up-to-date stacking code.
-In general, uppercase data arrays refer to all halos, and lower case data arrays refer to a specific halo
######## 
'''

## IMPORT MODULES ##
print '...importing modules'
import numpy as np
import astropy.io.fits as fits
from numpy.linalg import norm
import matplotlib.pyplot as mp
import astStats
import sys
import time
import cPickle as pkl
from mpl_toolkits.mplot3d import Axes3D

from caustic_class_stack2D import *
from caustic_universal_stack2D import *
from CausticMass import Caustic,CausticSurface,MassCalc


## FLAGS ##
self_stack	= False				# Run self-stack or bin-stack
scale_data	= True				# Scale data by r200 and vdisp if True
use_flux	= True				# Using Flux if True, using Sophie if False
write_data 	= True				# Write Data to Result directories if True
light_cone	= False				# Input RA|DEC projection data if True, if False inputting x,y,z 3D data
clean_ens	= False				# Do an extra shiftgapper on ensemble before the lines of sight get stacked.
small_set	= False				# 100 Halo Set or 2000 Halo Set
run_los		= False				# Run caustic technique on each line of sight?

## CONSTANTS ##
c 		= 2.99792e5			# speed of light in km/s
h		= 1.0				# Hubble Constant, unitless
H0		= h*100.0			# Hubble Constant, km s-1 Mpc-1
q		= 50.0				# Scale of Gaussian Kernel Density Estimator
beta		= 0.2				# Velocity Anisotropy Beta parameter, if constant profile
fbeta		= 0.65				# fbeta value, see 'Diaferio 1999'
r_limit 	= 1.5				# Radius Cut Scaled by R200
v_limit		= 3500.0			# Velocity Cut in km/s
data_set	= 'Guo30_2'			# Data set to draw semi analytic data from
halo_num	= 2100				# Total number of halos loaded
run_time	= time.asctime()		# Time when program was started

## RUN DEPENDENT CONSTANTS ##
run_num		= int(sys.argv[1])		# Run Number: referring to the run_num th iteration of this program via job array in FLUX
clus_num	= int(sys.argv[2])		# Number of Ensembles to build and solve for in this run
gal_num		= int(sys.argv[3])		# Number of galaxies taken per line of sight
line_num	= int(sys.argv[4])		# Number of lines of sight to stack over
method_num	= int(sys.argv[5])		# Ensemble Build Method Number
cell_num	= sys.argv[6]			# Cell Number ID corresponding to given gal_num & line_num geometry in a Run Table
table_num	= int(sys.argv[7])		# Table Re-Run Version	
data_loc	= 'binstack_run_table'+str(table_num)	# Parent Directory where write_loc directories live
try: 
	run_los = bool(int(sys.argv[8]))	# If fed 8th arg value as True, run_los
except:
	pass

if use_flux == True: 
	root=str('/nfs/christoq_ls')		# Change directory scheme if using flux or sophie
else: 
	root=str('/n/Christoq1')

if self_stack == True:								# Change Write Directory Depending on Parameters
	write_loc = 'ss_m'+str(method_num)+'_run'+str(cell_num)			# Self Stack data-write location
	stack_range = np.arange(run_num*clus_num,run_num*clus_num+clus_num)	# Range of halos, each to be stacked individually
else:
	write_loc = 'bs_m'+str(method_num)+'_run'+str(cell_num)			# Bin Stack data-write location
	stack_range = np.arange(run_num*clus_num*line_num,run_num*clus_num*line_num+clus_num*line_num)


## Make dictionary for above constants
varib = {'c':c,'h':h,'H0':H0,'q':q,'beta':beta,'fbeta':fbeta,'r_limit':r_limit,'v_limit':v_limit,'data_set':data_set,'halo_num':halo_num,'gal_num':gal_num,'line_num':line_num,'method_num':method_num,'write_loc':write_loc,'data_loc':data_loc,'root':root,'self_stack':self_stack,'scale_data':scale_data,'use_flux':use_flux,'write_data':write_data,'light_cone':light_cone,'run_time':run_time,'clean_ens':clean_ens,'small_set':small_set,'run_los':run_los,'run_num':run_num,'clus_num':clus_num,'cell_num':cell_num,'stack_range':stack_range}

## INITIALIZATION ##
U = universal(varib)
C = Caustic()
S = Stack(varib,U,C,CausticSurface,MassCalc)		# Pass varib dictionary, and pointer to other classes.
SS = SelfStack(varib,U,C,S,CausticSurface,MassCalc)	# Pass varib dictionary, and pointer to other classes and Stack instance
BS = BinStack(varib,U,C,S,CausticSurface,MassCalc)	# Pass varib dictionary, and pointer to other classes and Stack instance.

###################
##### PROGRAM #####
###################
U.print_separation('## Running caustic_mass_stack2D.py')
U.print_varibs(varib)

## Load Halo Data
U.print_separation('# ...Loading Halos',type=2)
HaloID,HaloData = U.load_halos()
# Sort Halos by A Priori Known Descending Mass (Mass Critical 200)
HaloID,HaloData = U.sort_halos(HaloID,HaloData)
# Unpack HaloData array into local namespace
M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ = HaloData

## Load Galaxy Data
U.print_separation('# ...Loading Galaxies',type=2)
Halo_P,Halo_V,Gal_P,Gal_V,G_Mags,R_Mags,I_Mags,HaloData = U.configure_galaxies(HaloID,HaloData)

# Get Gal_P in physical coordinates
Gal_P2 = []
if self_stack == True:
	for [i,k] in zip(np.arange(clus_num),stack_range):
		Gal_P2.append((Gal_P[i].T-Halo_P[k]).T)
else:
	for [i,l] in zip(np.arange(clus_num*line_num),stack_range):
		Gal_P2.append((Gal_P[i].T-Halo_P[l]).T)

# If Bin Stacking:
#	- Create Ensemble R200 and HVD arrays
#	- Change order of halos to bin upon
#	- Create any other arrays needed
if self_stack == False:
	BinData = U.Bin_Calc(HaloData,varib)
	BIN_M200,BIN_R200,BIN_HVD = BinData

# Initialize Multi-Ensemble Array to hold resultant data
STACK_DATA = []

#  j: index for each ensemble, w.r.t. the FOR loop
#  k: index for each final ensemble, w.r.t. total number of final ensembles
#  l: index for line of sight, w.r.t. total lines of sight of run
# (s: index for each to-be-stacked halo w.r.t. HaloData arrays)
U.print_separation('# ...Starting Ensemble Loop',type=2)
for j in range(clus_num):	# iterate over # of ensembles to build and solve for

	# Total Halo Index
	k = run_num*clus_num + j

	# Build Ensemble and Run Caustic Technique
	if self_stack:
		stack_data = SS.self_stack_clusters(HaloID,HaloData,Halo_P,Halo_V,Gal_P,Gal_V,G_Mags,R_Mags,I_Mags,k,j)
	else:
		stack_data = BS.bin_stack_clusters(HaloID,HaloData,BinData,Halo_P,Halo_V,Gal_P,Gal_V,G_Mags,R_Mags,I_Mags,k,j)

	# Unpack data
	ens_r,ens_v,ens_gal_id,ens_gmags,ens_rmags,ens_imags,ens_hvd,ens_caumass,ens_caumass_est,ens_causurf,ens_nfwsurf,los_r,los_v,los_gal_id,los_gmags,los_rmags,los_imags,los_hvd,los_caumass,los_caumass_est,los_causurf,los_nfwsurf,x_range,sample_size,pro_pos = stack_data

	if self_stack == True:
		# Get 3D data
		gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d = U.get_3d(Gal_P2[j],Gal_V[j],G_Mags[j],R_Mags[j],I_Mags[j],ens_gmags,ens_rmags,ens_imags)		
	else:	# I don't yet know how to do this efficiently for bin stacking
		gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d = [],[],[],[],[],[]

	# Combine into stack_data
	stack_data = [ens_r,ens_v,ens_gal_id,ens_gmags,ens_rmags,ens_imags,ens_hvd,ens_caumass,ens_caumass_est,ens_causurf,ens_nfwsurf,los_r,los_v,los_gal_id,los_gmags,los_rmags,los_imags,los_hvd,los_caumass,los_caumass_est,los_causurf,los_nfwsurf,x_range,sample_size,pro_pos,gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d]

	# Append to STACK_DATA
	STACK_DATA.append(stack_data)


# Finished Loop
U.print_separation('#...Finished Ensemble Loop',type=2)

	
### Save Data into Fits Files ###
if write_data == True:
	U.print_separation('#...Starting Data Write',type=2)
	for m in range(clus_num):
		n = run_num*clus_num + m
		U.print_separation("Writing Data for Ensemble #"+str(n),type=2)
		pkl_file = open(root+'/nkern/Stacking/'+data_loc+'/'+write_loc+'/Ensemble_'+str(n)+'_Data.pkl','wb')
		output = pkl.Pickler(pkl_file)
		output.dump(STACK_DATA[m])
		output.dump(varib)
		pkl_file.close()
	
	U.print_separation('#...Finished Data Write',type=2)


## Wait until at least 60 seconds is up
duration = (float(time.asctime()[11:13])*3600+float(time.asctime()[14:16])*60+float(time.asctime()[17:19])) - (float(run_time[11:13])*3600+float(run_time[14:16])*60+float(run_time[17:19]))
if duration < 60:
	time.sleep(60-duration)
	

U.print_separation('## Finished caustic_mass_stack2D.py'+'\n'+'Start:'+'\t'+run_time+'\n'+'End:'+'\t'+time.asctime(),type=1)


