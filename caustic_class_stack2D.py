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
from CausticMass import *
import numpy.random as npr
import astStats

## Program ##

class selfstack:

	def __init__(self,varib):
		''' Initial function for class selfstack '''
		# Adding dictionary varib to class namespace
		self.__dict__.update(varib)
		# Initializing the universal class 
		self.U = universal(varib)			
		self.C = Caustic()
		self.CS = CausticSurface()
		self.MC = MassCalc()
	
	def build_ensemble(self,r,v,mags,halodata,l):
		''' 
		- This function selects specific galaxies per line of sight using a sepcified method of stacking
		- The current method of interloper treatment is using CausticMass.py's ShiftGapper technique
		- Certain Measures were taken when using the ShiftGapper to ensure gal_num galaxies within r200
		- Interloper treatment is always done for the LOS, and can be done for ensemble_los if desired
		'''
		##
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
			end = within[:gal_num + excess + 1][-1]		# instead of indexing I am slicing b/c of chance of not enough gals existing...	
			# Build Ensemble (en => ensemble)
			if clean_en_los == True:
				excess = gal_num * 2.0 / 5.0		# make excess a bit larger than previously defined
				end = within[:gal_num + excess + 1][-1]
				r2,v2,mags2 = self.C.shiftgapper(np.vstack([r[:end],v[:end],mags[:end]]).T).T # Shiftgapper inputs and outputs data as transpose...
				within = np.where(r2<r_crit200)[0]	# re-calculate within array with new sample
				excess = gal_num / 5.0
				end = within[:gal_num + excess + 1][-1]
				# Append to ensemble array
				en_r,en_v,en_m = r2[:end],v2[:end],mags2[:end]
			else:
				en_r,en_v,en_m = r[0:end],v[0:end],mags[0:end]
			# Build Line of Sight (ln => line of sight)
			# shiftgapper on line of sight
			excess = gal_num * 2.0 / 5.0
			end = within[:gal_num + excess + 1][-1]
			r2,v2,mags2 = self.C.shiftgapper(np.vstack([r[:end],v[:end],mags[:end]]).T).T
			within = np.where(r2<r_crit200)[0]		# re-calculate within array with new sample
			# Now feed ln arrays correct gal_num richness within r200
			end = within[:gal_num + 1][-1]
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
			excess = gal_num * 2.0 / 5.0
			samp_num = gal_num + excess			# sample size of randomly generated numbers, start too low on purpose, then raise in loop
			loop = True					# break condition
			while loop == True:				# see method 0 comments on variables such as 'excess' and 'within' and 'end'
				for j in range(3):			# before upping sample size, try to get a good sample a few times
					rando = npr.randint(0,samp_size,samp_num)
					within = np.where(r[rando]<=r_crit200)[0]
					if len(within) >= gal_num + excess:
						loop = False
				if len(within) < gal_num + excess:
					samp_num += 2
			### Build Ensemble
			if clean_en_los == True:
				r2,v2,mags2 = self.C.shiftgapper(np.vstack([r[rando],v[rando],mags[rando]]).T).T
				within = np.where(r2<r_crit200)[0]
				excess = gal_num / 5.0
				end = within[:gal_num + excess + 1][-1]
				# Append to ensemble array
				en_r,en_v,en_m = r2[:end],v2[:end],mags2[:end]
			else:
				excess = gal_num / 5.0
				end = within[:gal_num + excess + 1][-1]
				en_r,en_v,en_m = r[rando][:end],v[rando][:end],mags[rando][:end]

			### Build LOS
			excess = gal_num / 5.0
			end = within[:gal_num + excess + 1][-1]
			r2,v2,mags2 = self.C.shiftgapper(np.vstack([r[rando][:end],v[rando][:end],mags[rando][:end]]).T).T
			within = np.where(r2<r_crit200)[0]
			end = within[:gal_num + 1][-1]
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


	def kernel_caustic_masscalc(self,r,v,halodata,derived_hvd,k,l=None):
		'''
		- This function runs the basic procedure of the Caustic Technique from CausticMass.py by using:
		- Caustic.gaussian_kernel()
		- CausticSurface.main()
		- MassCalc.main()
		'''
		## Print Details
		if l == None:
			text = '## Working On Cluster #'+str(k)+''
		else:
			text = '## Working On Cluster #'+str(k)+'\n## Line of Sight #'+str(l)+''
		self.U.print_separation(text,type=2)


		## Unpack HaloData Array into Namespace
		m_crit200,r_crit200,z,srad,esrad,hvd = halodata

		## Mirror Phase Space Velocity Data if mirror == True
		mirror = True
		if mirror == True:
			rvalues,vvalues = np.append(r,r),np.append(v,-v)		
		else:
			rvalues,vvalues = r,v

		## Perform Kernel Density Estimation
		# Set Kernel Density Phase Space limits in Mpc and km/s
		xmax = 6.0			# Mpc
		ymax = 5000.0			# km/s
		self.C.gaussian_kernel(rvalues,vvalues,r_crit200,normalization=self.H0,scale=self.q,xmax=xmax,ymax=ymax,xres=200,yres=220)
		self.C.img_tot = self.C.img/np.max(np.abs(self.C.img))
 		self.C.img_grad_tot = self.C.img_grad/np.max(np.abs(self.C.img_grad))
 		self.C.img_inf_tot = self.C.img_inf/np.max(np.abs(self.C.img_inf))
		
		## Define Beta Profile
		beta_profile = False		# beta profile as a funct of radius?
		if beta_profile == True:
			xbeta,abeta = np.loadtxt(''+str(root)+'/nkern/Caustic/average_betaprofile.tab',dtype='float',usecols=(0,1),unpack=True)
			fit = np.polyfit((xbeta*r_crit200)[xbeta<4],abeta[xbeta<4],6)
			self.beta = fit[0]*self.C.x_range**6+fit[1]*self.C.x_range**5+fit[2]*self.C.x_range**4+fit[3]*self.C.x_range**3+fit[4]*self.C.x_range**2+fit[5]*self.C.x_range+fit[6]	
		else:				# constant beta w/ radius
			self.beta = np.zeros(self.C.x_range.size) + self.beta

		## Caustic Surface Calculation
		# This function takes RA, DEC and Z as first input, for now it is empty, all relavent files go to class dictionary
		self.CS.main(np.zeros(0),self.C.x_range,self.C.y_range,self.C.img_tot,r200=r_crit200,halo_scale_radius=srad,halo_scale_radius_e=esrad,halo_vdisp=derived_hvd,beta=self.beta)

		## Mass Calculation, leave clus_z blank for now
		self.MC.main(self.C.x_range,self.CS.vesc_fit,derived_hvd,0,r200=r_crit200,fbr=self.fbeta,H0=self.H0)

		return self.MC.M200,self.CS.Ar_finalD,self.CS.vesc_fit	



	def self_stack_clusters(self,HaloID,HaloData,Halo_P,Halo_V,Gal_P,Gal_V,Gal_Mags,k):
		''' Building Ensemble Cluster and Calculating Property Statistics '''
		## Unpack HaloData array 
		M_crit200,R_crit200,Z,SRAD,ESRAD,HVD = HaloData

		## Define Arrays for Building Ensemble and LOS
		# Ensemble Arrays:	[Successive Ensemble Number][Data]
		# Line of Sight Arrays:	[Line Of Sight][Data]
		ens_r, ens_v, ens_m, ens_hvd = [], [], [], []
		ens_caumass, ens_causurf, ens_nfwsurf = [], [], []
		los_r, los_v, los_m, los_hvd = [], [], [], []
		los_caumass, los_causurf, los_nfwsurf = [], [], []

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
				ln_hvd = np.std( np.copy(ln_v)[ln_within] )
			else:
				# Calculate hvd with astStats biweightScale (see Beers 1990)
				ln_hvd = astStats.biweightScale(np.copy(ln_v)[ln_within],9.0)

			# Run Caustic Technique for LOS mass estimation
			ln_caumass,ln_causurf,ln_nfwsurf = self.kernel_caustic_masscalc(r,v,HaloData.T[k],ln_hvd,k,l)
			
			# Append LOS Data Arrays
			los_r.append(ln_r)
			los_v.append(ln_v)
			los_m.append(ln_m)
			los_hvd.append(ln_hvd)
			los_caumass.append(ln_caumass)
			los_causurf.append(ln_causurf)
			los_nfwsurf.append(ln_nfwsurf)

		# Shiftgapper for Ensemble Interloper treatment
		ens_r,ens_v,ens_m = self.C.shiftgapper(np.vstack([ens_r,ens_v,ens_m]).T).T

		# Reduce system to gal_num richness within r200
		within = np.where(ens_r <= R_crit200[k])[0]
		end = within[:self.gal_num*self.line_num + 1][-1]
		ens_r = ens_r[:end]
		ens_v = ens_v[:end]
		ens_m = ens_m[:end]

		# Calculate HVD
		en_hvd = astStats.biweightScale(np.copy(ens_v),9.0)

		# Caustic Technique for Ensemble
		en_caumass,en_causurf,en_nfwsurf = self.kernel_caustic_masscalc(ens_r,ens_v,HaloData.T[k],en_hvd,k)

		# Append Ensemble Data Arrays
		ens_hvd.append(en_hvd)
		ens_caumass.append(en_caumass)

		# Turn into numpy arrays
		ens_r,ens_v,ens_m = np.array(ens_r),np.array(ens_v),np.array(ens_m)
		ens_hvd,ens_caumass = np.array(ens_hvd),np.array(ens_caumass)
		ens_causurf,ens_nfwsurf = np.array(en_causurf),np.array(en_nfwsurf)
		los_r,los_v,losm = np.array(los_r),np.array(los_v),np.array(los_m)
		los_hvd,los_caumass = np.array(los_hvd),np.array(los_caumass)
		los_causurf,los_nfwsurf = np.array(los_causurf),np.array(los_nfwsurf)

		return ens_r,ens_v,ens_m,ens_hvd,ens_caumass,ens_causurf,ens_nfwsurf,los_r,los_v,los_m,los_hvd,los_caumass,los_causurf,los_nfwsurf,self.C.x_range	
	





