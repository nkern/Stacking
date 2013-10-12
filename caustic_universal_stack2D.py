#caustic_universal_stack2D.py
'''
This script contains the class universal, which is used by caustic_class_stack2D.py and caustic_mass_stack2D.py

universal:
	functions:
	attributes:
'''

## Import Modules ##
import CausticMass as cm
import numpy as np
import pyfits
import random
from scipy import weave
from scipy.weave import converters
import cosmolopy.distance as cd
import time

## Program ##

class universal:

	def __init__(self,varib):
		''' Adding permanent program variables to class namespace'''
		self.__dict__.update(varib)

	def load_halos(self):
		'''This function loads halo data and makes cosmology corrections'''
		HaloID = np.loadtxt(self.root+'/nkern/Caustic/biglosclusters.csv', delimiter=',', dtype='string', usecols=(0,), unpack=True)
		HPX,HPY,HPZ,HVX,HVY,HVZ = np.loadtxt(self.root+'/nkern/Caustic/biglosclusters.csv', delimiter=',', dtype='float', usecols=(9,10,11,12,13,14), unpack=True)
		SRAD,ESRAD,R_crit200,M_crit200,HVD,Z = np.loadtxt(self.root+'/nkern/Caustic/Millbig_concentrations.phys_phys.csv', delimiter=',', dtype='float', usecols=(1,2,5,7,9,12), unpack=True)
		# Hubble Constant Coefficient
		R_crit200,M_crit200,HPX,HPY,HPZ = R_crit200/self.h,M_crit200/self.h,HPX/self.h,HPY/self.h,HPZ/self.h
		# Cosmological Correction
		for l in xrange(len(HaloID)):	
			HPX[l],HPY[l],HPZ[l] = HPX[l]/(1+Z[l]),HPY[l]/(1+Z[l]),HPZ[l]/(1+Z[l])	
			# Fix weird SRAD values, if R_crit200/SRAD = Conc > 2, set SRAD=R_crit200/2
			if R_crit200[l]/SRAD[l] < 2.0:
				SRAD[l] = R_crit200[l] / 2.0
		#Ordering of halos in arrays is identical to biglosclusters' inherent ordering.
			
		return HaloID, np.vstack([M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ])


	def sort_halos(self,HaloID,HaloData):
		''' Sort Halo Data by some Criteria '''
		# Unpack Array HaloData into local namespace for easier use and clarity
		M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ = HaloData	
		# Sort Arrays by descending M_crit200	
		sort = np.argsort(M_crit200)[::-1]	
		HaloID = HaloID[sort]
		M_crit200 = M_crit200[sort]
		R_crit200 = R_crit200[sort]
		Z = Z[sort]
		SRAD = SRAD[sort]
		ESRAD = ESRAD[sort]
		HVD = HVD[sort]
		HPX = HPX[sort]
		HPY = HPY[sort]
		HPZ = HPZ[sort]
		HVX = HVX[sort]
		HVY = HVY[sort]
		HVZ = HVZ[sort]
		# Return packed array
		return HaloID, np.vstack([M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ]) 

	def configure_galaxies(self,HaloID,HaloData):
		''' Loads galaxy data from halo list, and converts to physical coordinates and corrects cosmological factors '''
		# Unpack Array HaloData into local namespace for easier use and clarity
		M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ = HaloData	

		Gal_Mags = []
		Gal_V = []
		Gal_P = []
		Halo_P = []
		Halo_V = []
		for k in range(self.halo_num):
			galdata,hpx,hpy,hpz,hvx,hvy,hvz = self.load_galaxies(HaloID[k],HaloData.T[k])
			# unpack array galdata into namespace
			gpx,gpy,gpz,gvx,gvy,gvz,gal_mags = galdata	
			halo_p	= np.array([hpx,hpy,hpz],float) 
			halo_v	= np.array([hvx,hvy,hvz],float)
			gal_p	= np.array([gpx,gpy,gpz],float)
			gal_v	= np.array([gvx,gvy,gvz],float)
			
			Gal_Mags.append(gal_mags)
			Halo_P.append(halo_p)
			Halo_V.append(halo_v)			
			Gal_P.append(gal_p)
			Gal_V.append(gal_v)			

		Halo_P,Halo_V = np.array(Halo_P),np.array(Halo_V)
		Gal_Mags = np.array(Gal_Mags)		

		return Halo_P,Halo_V,Gal_P,Gal_V,Gal_Mags,HaloData[0:6]

	def load_galaxies(self,haloid,halodata):
		''' Loads haloid galaxies from a local directory '''
		# Unpack array halodata into local namespace
		m_crit200,r_crit200,z,srad,esrad,hvd,hpx,hpy,hpz,hvx,hvy,hvz = halodata
		# Open semi analytics	
		f = pyfits.open(self.root+'/giffordw/Millenium/30Mpchalos/'+haloid+'.'+self.data_set+'.fits')
		data = f[1].data
		gal_z,gpx,gpy,gpz,gvx,gvy,gvz,mags = data.field(13),data.field(17),data.field(18),data.field(19),data.field(20),data.field(21),data.field(22),data.field(63)
		# Cosmology corrections
		gpx,gpy,gpz = (gpx/(1+z)/self.h),(gpy/(1+z)/self.h),(gpz/(1+z)/self.h)
		# convert to physical coordinates
		gvx,gvy,gvz = gvx-hvx,gvy-hvy,gvz-hvz
		mags = np.array(mags,float)
		# remove BCG from sample
		BCG = np.where(gpx != hpx)

		return np.vstack([ gpx[BCG], gpy[BCG], gpz[BCG], gvx[BCG], gvy[BCG], gvz[BCG], mags[BCG] ]),hpx,hpy,hpz,hvx,hvy,hvz


	def scale_gals(self,r,v,r_crit200,hvd):
		''' Scales galaxy projected radius and velocity data by a given radius and velocity dispersion '''
		r /= r_crit200
		v /= hvd
		return r, v

	def pick_pos(self,halo_p):
	        '''Picks a random position for the observer a given distance away from the center'''
		x = random.uniform(-1,1)
		y = random.uniform(-1,1)
		z = random.uniform(-1,1)
		unit = np.array([x,y,z])/(x**2+y**2+z**2)**(.5)
		# move the position randomly 30Mpc away
	        return halo_p + 30*unit

	def line_of_sight(self,gal_p,gal_v,halo_p,halo_v):
		'''Line of Sight Calculations to mock projected data, if given 3D data'''
		# Pick Position
		new_pos = self.pick_pos(halo_p) 

		# New Halo Information
		halo_dist = ((halo_p[0]-new_pos[0])**2 + (halo_p[1]-new_pos[1])**2 + (halo_p[2]-new_pos[2])**2)**0.5
		halo_pos_unit = np.array([halo_p[0]-new_pos[0],halo_p[1]-new_pos[1],halo_p[2]-new_pos[2]]) / halo_dist
		halo_vlos = np.dot(halo_pos_unit, halo_v)

		# New Galaxy Information
		gal_p = np.array(gal_p)
		gal_v = np.array(gal_v)
		gal_dist = ((gal_p[0]-new_pos[0])**2 + (gal_p[1]-new_pos[1])**2 + (gal_p[2]-new_pos[2])**2)**0.5
		gal_vlos = np.zeros(gal_dist.size)
		gal_pos_unit = np.zeros((3,gal_dist.size))	#vector from new_p to gal	
		n = gal_dist.size
		# Line of sight
		code = """
		int u,w;
		for (u=0;u<n;++u){
		for(w=0;w<3;++w){
		gal_pos_unit(w,u) = (gal_p(w,u)-new_pos(w))/gal_dist(u);
		}
		gal_vlos(u) = gal_pos_unit(0,u)*gal_v(0,u)+gal_pos_unit(1,u)*gal_v(1,u)+gal_pos_unit(2,u)*gal_v(2,u);
		}
		"""
		fast = weave.inline(code,['gal_pos_unit','n','gal_dist','gal_vlos','gal_v','new_pos','gal_p'],type_converters=converters.blitz,compiler='gcc')
		angles = np.arccos(np.dot(halo_pos_unit,gal_pos_unit))
		r = angles*halo_dist
		v_pec = gal_vlos-halo_vlos*np.dot(halo_pos_unit,gal_pos_unit)
		z_clus_cos = self.H0*halo_dist/self.c
		z_clus_pec = 0#halo_vlos/c
		z_clus_obs = (1+z_clus_pec)*(1+z_clus_cos)-1
		z_gal_cos = self.H0*gal_dist/self.c
		z_gal_pec = gal_vlos/self.c
		z_gal_obs = (1+z_gal_pec)*(1+z_gal_cos)-1
		v = self.c*(z_gal_obs-z_clus_obs)/(1+z_clus_obs)
		#gal_vdisp3d[i] = np.sqrt(astStats.biweightScale(gal_v[0][np.where(gal_radius<=HaloR200[i])]-Halo_V[0],9.0)**2+astStats.biweightScale(gal_v[1][np.where(gal_radius<=HaloR200[i])]-Halo_V[1],9.0)**2+astStats.biweightScale(gal_v[2][np.where(gal_radius<=HaloR200[i])]-Halo_V[2],9.0)**2)/np.sqrt(3)
		#print 'MY VELOCITY OF GALAXIES', gal_vdisp3d[i]
#		particle_vdisp3d[i] = HVD*np.sqrt(3)
#		gal_rmag_new = gal_abs_rmag# + 5*np.log10(gal_dist*1e6/10.0)

		return r, v

	def limit_gals(self,r,v,mags,r200,hvd):
		''' Sort data by magnitude, and elimite values outside phase space limits '''
		# Sort by ascending magnitude (bright to dim)
		sorts = np.argsort(mags)
		r,v,mags = r[sorts],v[sorts],mags[sorts]
		# Limit Phase Space
		sample = np.where( (r < r200*self.r_limit) & (v > -self.v_limit) & (v < self.v_limit) )[0] 
		r,v,mags = r[sample],v[sample],mags[sample]
	
		return r,v,mags

	def print_varibs(self,varib):
		print "Start Time		=",time.asctime()
		print "halo_num		=",varib['halo_num']
		print "ens_num			=",varib['ens_num']
		print "gal_num			=",varib['gal_num']
		print "line_num		=",varib['line_num']
		print "method_num		=",varib['method_num']
		print "write_loc		=",varib['write_loc']
		print "data_set		=",varib['data_set']
		print "self_stack		=",varib['self_stack']
		print "write_data		=",varib['write_data']
		print "scale_data		=",varib['scale_data']
		print "light_cone		=",varib['light_cone']
		return	

	def print_separation(self,text,type=1):
		if type==1:
			print ''
			print '-'*60
			print str(text)
			print '-'*60
		elif type==2:
			print ''
			print str(text)
			print '-'*30
		return




