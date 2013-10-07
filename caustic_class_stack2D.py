## caustic_class_stack2D.py
'''
#########
- This script contains the class 'selfstack', used by caustic_mass_stack2D.py.
- Parts of the program are either modelled or directly taken from D. Gifford's CausticMass.py code. 
- This is the most up-to-date stacking code.
#########

self_stack:
	functions:
	attributes:
'''

## Import Modules ##
import numpy as np
from caustic_universal_stack2D import universal
from CausticMass import Caustic 
import numpy.random as npr
import astStats

## Program ##

class selfstack(object):

	def __init__(self,varib):
		''' Initial function for class selfstack '''
		# Adding dictionary varib to class namespace
		self.__dict__.update(varib)
		# Initializing the universal class 
		self.U = universal(varib)			
		self.C = Caustic()
	
	def build_ensemble(self,r,v,mags,halodata,l):
		''' 
		- This function selects specifc galaxies per line of sight using a sepcified method of stacking
		- The current method of interloper treatment is using CausticMass.py's ShiftGapper technique
		- Certain Measures were taken when using the ShiftGapper to ensure gal_num galaxies within r200
		- anterloper treatment is always done for the LOS, and can be done for ensemble_los if desired
		'''
		## Constants and Flags
		clean_en_los = False					# This flag toggles extra shiftgapper before ensemble stacking
		gal_num = self.gal_num	
	
		## Program
		# Unpack halodata array into local namespace
		m_crit200,r_crit200,z,srad,esrad,hvd = halodata	
		# Sort galaxies by Magnitude
		bright = np.argsort(mags)
		r,v,mags = r[bright],v[bright],mags[bright]

		if self.method_num == 0:
			'''Picking top brightest galaxies, such that there are gal_num galaxies within r200'''
			# define indicies of galaxies within r200
			within = np.where(r<r_crit200)[0]
			# pick out gal_num 'th index in list, (include extra to counter shiftgapper's possible dimishement of richness)
			excess = gal_num / 5.0				# galaxy excess to counter shiftgapper
			end = within[gal_num + excess]		
			# Build Ensemble (en => ensemble)
			if clean_en_los == True:
				excess = gal_num * 2.0 / 5.0		# make excess a bit larger than previously defined
				end = within[gal_num + excess]
				r2,v2,mags2 = self.C.shiftgapper(np.vstack([r[:end],v[:end],mags[:end]]).T).T # Shiftgapper inputs and outputs data as transpose...
				within = np.where(r2<r_crit200)[0]	# re-calculate within array with new sample
				excess = gal_num / 5.0
				end = within[gal_num + excess]
				# Append to ensemble array
				en_r,en_v,en_m = r2[:end],v2[:end],mags2[:end]
			else:
				en_r,en_v,en_m = r[0:end],v[0:end],mags[0:end]
			# Build Line of Sight (ln => line of sight)
			# shiftgapper on line of sight
			excess = gal_num * 2.0 / 5.0
			end = within[gal_num + excess]
			r2,v2,mags2 = self.C.shiftgapper(np.vstack([r[:end],v[:end],mags[:end]]).T).T
			within = np.where(r2<r_crit200)[0]		# re-calculate within array with new sample
			# Now feed ln arrays correct gal_num richness within r200
			end = within[gal_num]
			ln_r,ln_v,ln_m = r2[:end],v2[:end],mags2[:end]	
			# Done! Now we have en_r and ln_r arrays, which will either be stacked (former) or put straight into Caustic technique (latter)
			# Note that the pure los and ensemble los do not necessarily have identical phase space due to shiftgapper treatment

		elif self.method_num == 1:
			'''Randomly choosing galaxies until gal_num galaxies are within r200'''
			# reduce size of sample to something reasonable within magnitude limits
			sample = gal_num * 25				# arbitrary coefficient, see sites page post Apr 24th, 2013 for more info
			r,v,mags = r[:sample],v[:sample],mags[:sample]
			samp_size = len(r)				# actual size of sample (might be less than gal_num*25)
			# create random numbered array for galaxy selection
			while True:					# see method 0 comments on variables such as 'excess' and 'within' and 'end'
				excess = gal_num * 2.0 / 5.0
				end = gal_num * 1.5
				rando = npr.randint(0,samp_size,end)
				within = np.where(r[rando]<=r_crit200)[0]
				if len(within) < gal_num + excess:
					end += 10
				else:
					break
			### Build Ensemble
			if clean_en_los == True:
				r2,v2,mags2 = self.C.shiftgapper(np.vstack([r[rando],v[rando],mags[rando]]).T).T
				within = np.where(r2<r_crit200)[0]
				excess = gal_num / 5.0
				end = within[gal_num + excess]
				# Append to ensemble array
				en_r,en_v,en_m = r2[:end],v2[:end],mags2[:end]
			else:
				excess = gal_num / 5.0
				end = within[gal_num + excess]
				en_r,en_v,en_m = r[rando][:end],v[rando][:end],mags[rando][:end]

			### Build LOS
			excess = gal_num / 5.0
			end = within[gal_num + excess]
			r2,v2,mags2 = self.C.shiftgapper(np.vstack([r[rando][:end],v[rando][:end],mags[rando][:end]]).T).T
			within = np.where(r2<r_crit200)[0]
			end = within[gal_num]
			ln_r,ln_v,ln_m = r2[:end],v2[:end],mags2[:end]
			# Done! Now we have en_r and ln_r arrays (ensemble and line of sight arrays)

		elif self.method_num == 2:
			'''Ordered Set of galaxies with respect to magnitude'''
			'''
			### NOT FINSHED YET !! ###
			### Build Ensemble 
			within = np.where(r<r_crit200)[0]
			excess = gal_num / 5.0
			end = gal_num + excess
			if clean_end_los == True:
				pass		
	
			end = within[gal_num*line_num + line_num]
			en_r,en_v,en_mags,en_gpx3d,en_gpy3d,en_gpz3d,en_gvx3d,en_gvy3d,en_gvz3d = r[:end][l::line_num],v[:end][l::line_num],mags[:end][l::line_num],gpx3d[:end][l::line_num],gpy3d[:end][l::line_num],gpz3d[:end][l::line_num],gvx3d[:end][l::line_num],gvy3d[:end][l::line_num],gvz3d[:end][l::line_num]	
			### Build LOS
			end = within[:gal_num*line_num + 5*line_num][-1] + 1
			rt,vt,magst,gpx3dt,gpy3dt,gpz3dt,gvx3dt,gvy3dt,gvz3dt = U.shiftgapper(vstack((r[:end][l::line_num],v[:end][l::line_num],mags[:end][l::line_num],gpx3d[:end][l::line_num],gpy3d[:end][l::line_num],gpz3d[:end][l::line_num],gvx3d[:end][l::line_num],gvy3d[:end][l::line_num],gvz3d[:end][l::line_num])).T).T
			within = where(rt<r_crit200)[0]
			end = within[:gal_num][-1] + 1
			r,v,mags,gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d = rt[:end],vt[:end],magst[:end],gpx3dt[:end],gpy3dt[:end],gpz3dt[:end],gvx3dt[:end],gvy3dt[:end],gvz3dt[:end] 
			###########################
			'''
			pass

		return en_r,en_v,en_m,ln_r,ln_v,ln_m


	def kernel_caustic_masscalc(self):
		'''
		- This function runs the basic procedure of the Caustic Technique from CausticMass.py by using:
		- Caustic.gaussian_kernel()
		- CausticSurface.main()
		- MassCalc.main()
		'''
		pass
		





	def self_stack_clusters(self,HaloID,HaloData,Halo_P,Halo_V,Gal_P,Gal_V,Gal_Mags,k):
		''' Building Ensemble Cluster and Calculating Property Statistics '''
		## Unpack HaloData array 
		M_crit200,R_crit200,Z,SRAD,ESRAD,HVD = HaloData

		## Define Arrays for Building Ensemble and LOS
		# Ensemble Arrays:	[Successive Ensemble Number][Data]
		# Line of Sight Arrays:	[Line Of Sight][Data]
		ens_r, ens_v, ens_m = [], [], []
		los_r, los_v, los_m = [], [], []
		## Loop over lines of sight
		for l in range(self.line_num):
			if self.light_cone == True:
				# Configure RA, DEC and Z into cluster-centric radius and velocity
				pass
			else:
				# Line of Sight Calculation for naturally 3D data
				r, v = self.U.line_of_sight(Gal_P[k],Gal_V[k],Halo_P[k],Halo_V[k])

			# Limit Data in Phase Space
			r, v, mags = self.U.limit_gals(r,v,Gal_Mags[k],R_crit200[k],HVD[k])					
			
			# Build LOS and Ensemble, with given method of stacking
			en_r,en_v,en_m,ln_r,ln_v,ln_m = self.build_ensemble(r,v,mags,HaloData.T[k],l)	

			# Build Ensemble Arrays
			ens_r.extend(en_r)
			ens_v.extend(en_v)
			ens_m.extend(en_m)
		
			# Calculate LOS HVD (this is after shiftgapper)
			ln_within = np.where(ln_r<R_crit200[k])[0]
			gal_count = len(ln_within)
			if gal_count <= 3:
				'''biweightScale can't take less than 4 elements'''
				# Calculate hvd with numpy std of galaxies within r200 (b/c this is quoted richness)
				los_hvd = np.std( copy(ln_v)[ln_within] )
			else:
				# Calculate hvd with astStats biweightScale (see Beers 1990)
				los_hvd = astStats.biweightScale(np.copy(ln_v)[ln_within],9.0)

			# Running LOS mass estimation
		#	self.kernel_caustic_masscalc(
			


	





