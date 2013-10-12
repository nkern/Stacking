## caustic.mass_stack2D.py
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
import pyfits
from numpy.linalg import norm
import matplotlib.pyplot as mp
import astStats
import sys
import time
import cPickle as pkl

from caustic_class_stack2D import *
from caustic_universal_stack2D import *
from CausticMass import Caustic,CausticSurface,MassCalc


## FLAGS ##
self_stack	= True			# Run self-stack or bin-stack
scale_data	= False			# Scale data by r200 and vdisp if True
use_flux	= True			# Using Flux if True, using Sophie if False
write_data 	= True			# Write Data to Result directories if True
light_cone	= False			# Input RA|DEC projection data if True, if False inputting x,y,z 3D data
one_ens		= True			# Only solve for one ensemble cluster if true, this is generally the case when using an HPC

## CONSTANTS ##
c 		= 2.99792e5		# speed of light in km/s
h		= 1.0			# Hubble Constant, unitless
H0		= h*100.0		# Hubble Constant, km s-1 Mpc-1
q		= 10.0			# Scale of Gaussian Kernel Density Estimator
beta		= 0.2			# Velocity Anisotropy Beta parameter, if constant profile
fbeta		= 0.65			# fbeta value, see 'Diaferio 1999'
r_limit 	= 1.5			# Radius Cut Scaled by R200
v_limit		= 3500.0		# Velocity Cut in km/s
data_set	= 'Guo30_2'		# Data set to draw semi analytic data from
halo_num	= 100			# Total number of halos loaded
run_time	= time.asctime()	# Time when program was started

## RUN DEPENDENT CONSTANTS ##
if len(sys.argv) > 1:				# If you feed the run with parameters
	ens_num		= int(sys.argv[1])	# Number of Ensembles to solve for
	gal_num		= int(sys.argv[2])	# Number of galaxies taken per line of sight
	line_num	= int(sys.argv[3])	# Number of lines of sight to stack over
	run_num		= int(sys.argv[4])	# Run Number ID corresponding to given gal_num & line_num geometry
	method_num	= int(sys.argv[5])	# Self Stacking Method to Use
else:
	ens_num		= 1
	gal_num		= 10
	line_num	= 10
	run_num		= 1
	method_num	= 1

if use_flux == True: 
	root=str('/nfs/christoq_ls')	# Change directory scheme if using flux or sophie
else: 
	root=str('/n/Christoq1')

if self_stack == True:							# Change Write Directory Depending on Parameters
	write_loc = 'ss_m'+str(method_num)+'_run'+str(run_num)		# Self Stack data-write location
else:
	write_loc = 'bs_m'+str(method_num)+'_run'+str(run_num)		# Bin Stack data-write location

## Make dictionary for above constants
varib = {'c':c,'h':h,'H0':H0,'q':q,'beta':beta,'fbeta':fbeta,'r_limit':r_limit,'v_limit':v_limit,'data_set':data_set,'halo_num':halo_num,'ens_num':ens_num,'gal_num':gal_num,'line_num':line_num,'method_num':method_num,'write_loc':write_loc,'root':root,'self_stack':self_stack,'scale_data':scale_data,'use_flux':use_flux,'write_data':write_data,'light_cone':light_cone,'one_ens':one_ens,'run_time':run_time}

## INITIALIZATION ##
U = universal(varib)
C = Caustic()
CS = CausticSurface()
MC = MassCalc()
SS = selfstack(varib,U,C,CS,MC)		# Pass selfstack class variable dictionary, and pointer to other class instances

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
Halo_P,Halo_V,Gal_P,Gal_V,Gal_Mags,HaloData = U.configure_galaxies(HaloID,HaloData)

## Solve for Ensemble number: ens_num 
U.print_separation('# ...Starting Ensemble Loop',type=2)
j = 0
for k in np.array([ens_num]):

	# Build Ensemble and Run Caustic Technique
	if self_stack:
		stack_data = SS.self_stack_clusters(HaloID,HaloData,Halo_P,Halo_V,Gal_P,Gal_V,Gal_Mags,k)
		# unpack data
		ens_r,ens_v,ens_m,ens_hvd,ens_caumass,ens_causurf,ens_nfwsurf,los_r,los_v,los_m,los_hvd,los_caumass,los_causurf,los_nfwsurf,x_range = stack_data

	else:
		BS.bin_stack_clusters()

	# Condense into Arrays
	

	j += 1




### Save Data into Pickle Files ###
if write_data == True:
	pkl_file = open(root+'/nkern/Stacking/stack_data/'+write_loc+'/Data.pkl','wb')
	output = pkl.Pickler(pkl_file)
	output.dump(varib)
	output.dump([HaloID,Halo_P,Halo_V,Gal_P,Gal_V,Gal_Mags,HaloData])
	output.dump(stack_data)










