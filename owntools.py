import numpy as np
import c2raytools as c2t
from scipy.signal import argrelextrema
from sklearn.cluster import KMeans
import glob

def threshold_21cm(cube, p=0.25):
	out = np.zeros(cube.shape)
	th  = (1.-p)*cube.min() + p*cube.max()
	out[cube<=th] = 1.
	return out, th

def threshold_constant_mean(cube):
	"""
	The input is the ionization fraction cube.
	"""
	xc = cube.mean()
	lin = np.squeeze(cube.reshape(-1,1))
	lin.sort()
	Nn = round((1-xc)*len(lin))
	xth = lin[Nn]
	array = np.zeros(cube.shape)
	array[cube>=xth] = 1
	print "The output contains a tuple with binary-cube and determined-threshold."
	return array, xth

def threshold_PDF_split(cube, bins=200, upper_lim=False):
	"""
	The input is the brightness temperature cube.
	"""
	array = np.zeros(cube.shape)
	ht = np.histogram(np.squeeze(cube.reshape(-1,1)), bins=bins)
	minm = argrelextrema(ht[0], np.less)
	t_th = (ht[1][minm[0][0]]+ht[1][minm[0][0]+1])/2.
	if upper_lim: array[cube<=t_th] = 1
	else: array[cube>=t_th] = 1
	print "The output contains a tuple with binary-cube and determined-threshold."
	return array, t_th

def threshold_kmeans(cube, upper_lim=False, mean_remove=True, n_jobs=1):
	"""
	The input is the brightness temperature cube.
	"""
	array = np.zeros(cube.shape)
	#km = KMeans(n_clusters=2)
	if mean_remove:
		if upper_lim: X = cube[cube<=cube.mean()].reshape(-1,1)
		else: X = cube[cube>=cube.mean()].reshape(-1,1)
	else:
	 	X  = cube.reshape(-1,1)
	y = KMeans(n_clusters=2, n_jobs=n_jobs).fit_predict(X)
	if X[y==1].mean()>X[y==0].mean(): t_th = X[y==0].max()
	else: t_th = X[y==1].max()
	if upper_lim: array[cube<=t_th] = 1
	else: array[cube>t_th] = 1
	print "The output contains a tuple with binary-cube and determined-threshold."
	return array, t_th
	
def threshold_kmeans_3cluster(cube, upper_lim=False, n_jobs=1):
	"""
	The input is the brightness temperature cube.
	"""
	km = KMeans(n_clusters=3, n_jobs=n_jobs)
	X  = cube.reshape(-1,1)
	array = np.zeros(X.shape)
	km.fit(X)
	y = km.labels_
	centers = km.cluster_centers_
	if upper_lim: true_label = centers.argmin()
	else: true_label = centers.argmax()
	array[y==true_label] = 1
	array = array.reshape(cube.shape)
	return array

def get_binary(array, thres, greater=False):
	binary = np.zeros(array.shape)
	if greater:
		binary[array <= thres] = 1 
	else:
		binary[array >= thres] = 1
	return binary 

	
def get_distribution(array, resolution=1.0, bins=100, sizes=False):
	if sizes:
		sizes = array
	else:
		mn, mx = array.min(), array.max()
		sizes  = np.arange(mx)+1.
		print mx
		for i in xrange(int(mx)):
			label = i+1
			sizes[i] = len(array[array==label])
			print label,
	ht   = np.histogram(np.log(sizes), bins=bins)
	vols, dist = np.zeros(len(ht[0])+1), np.zeros(len(ht[0])+1)
	vols      = np.exp(ht[1])*resolution
	dist[:-1] = ht[0]
	dist[dist==0] = dist[dist!=0].min()
	return sizes, dist/np.sum(dist), vols

def plot_fof_sizes(sizes, bins=100):
	lg = np.log10(np.array(sizes))
	ht = np.histogram(lg, bins=bins)
	xx = 10**ht[1]
	yy = ht[0]*xx[:-1]
	zz = yy/np.sum(yy)
	dummy = zz[zz!=0].min()/1000.
	zz[zz==0] = dummy
	zz = np.hstack((zz,dummy))
	print "The output is Size, Size**2 dP/d(Size), lowest value"
	return xx, zz, dummy

def draw_sphere(array, centre, radius, label=1, periodic=True):
	nx, ny, nz = array.shape
	rang = np.arange(-radius, radius+1)
	for i in rang:
		for j in rang:
			for k in rang:
				x = centre[0]+i
				y = centre[1]+j
				z = centre[2]+k
				range_condition = x in range(0,nx) and y in range(0,ny) and z in range(0,nz)
				if range_condition or periodic:
					if i*i+j*j+k*k <= radius*radius:
						if x>=nx: x = x-nx
						if y>=ny: y = y-ny
						if z>=nz: z = z-nz
						array[x,y,z] = label
	if periodic: print "Periodic sphere of radius", radius, "made at", centre
	else: print "Non-periodic sphere of radius", radius, "made at", centre

def draw_circle(array, centre, radius, label=1, periodic=True):
	nx, ny = array.shape
	rang = np.arange(-radius, radius+1)
	for i in rang:
		for j in rang:
			x = centre[0]+i
			y = centre[1]+j
			range_condition = x in range(0,nx) and y in range(0,ny)
			if range_condition or periodic:
				if i*i+j*j <= radius*radius: array[x,y] = label
	if periodic: print "Periodic circle of radius", radius, "made at", centre
	else: print "Non-periodic circle of radius", radius, "made at", centre
				
				
def evolving_sphere(R, z, res=256, c_res=c2t.conv.LB, alpha=4.):
	#assert(res==2*Rb)
	ar = np.zeros((res,res,res))
	cdist  = c2t.z_to_cdist(z)
	cdist_low  = cdist-c2t.conv.LB/2
	cdist_high = cdist+c2t.conv.LB/2
	cdists = np.linspace(cdist_low, cdist_high, res)
	zs = c2t.cdist_to_z(cdists)
	cx, cy, cz = np.round(res/2.), np.round(res/2.), np.round(res/2.)
	RR = R*np.exp(alpha*(z-zs))   # Exponential growth
	#RR = R+alpha*R*(z-zs)	       # Linear
	#RR = R+alpha*R*(z-zs)**3      # Cubic
	for k in xrange(res):
		r = RR[k]*(1+z)/(1+zs[k])
		x = np.abs(k-cz)
		if r>x:
			rp = np.sqrt(r**2-x**2)
			dum = np.zeros((res,res))
			draw_circle(dum, [cx,cy], rp, label=1, periodic=True)
			#dum[cx,cy+rp], dum[cx,cy-rp] = 1, 1
			ar[:,:,k] = dum
	return ar, RR[0]*(1+z)/(1+zs[-1]), RR[-1]*(1+z)/(1+zs[-1])
			
def coeval_21cm(xfrac_dir, dens_dir, z, interpolation='linear'):
	"""
	xfrac_dir     : Give the path that contains the xfrac-files.
	dens_dir      : Give the path that contains the density-files.
	z	      : Redshift.
	interpolation : This is used when the coveal cube at that redshift is not available.
	"""
	xfrac = coeval_xfrac(xfrac_dir, z, interpolation=interpolation)
	dens  = coeval_dens(dens_dir, z, interpolation=interpolation)
	return c2t.calc_dt(xfrac, dens, z=z)


def coeval_xfrac(xfrac_dir, z, interpolation='linear'):
	"""
	xfrac_dir     : Give the path that contains the xfrac-files.
	z	      : Redshift.
	interpolation : This is used when the coveal cube at that redshift is not available.
	"""
	if not interpolation in ['linear']: #, 'step', 'sigmoid', 'step_cell'
		raise ValueError('Unknown interpolation type: %s' % interpolation)
	xfrac_files = glob.glob(xfrac_dir + '/xfrac3d_*.bin')
	xfrac_zs = None
	xfrac_zs = c2t.lightcone._get_file_redshifts(xfrac_zs, xfrac_files)
	if z in xfrac_zs:
		xfrac = c2t.XfracFile(xfrac_files[np.argwhere(z==xfrac_zs)]).xi
	else:
		z_l = xfrac_zs[xfrac_zs<z].max()
		z_h = xfrac_zs[xfrac_zs>z].min()
		xfrac_l = c2t.XfracFile(xfrac_files[xfrac_zs[xfrac_zs<z].argmax()]).xi
		xfrac_h = c2t.XfracFile(xfrac_files[xfrac_zs[xfrac_zs>z].argmin()]).xi
		xfrac = xfrac_h + (xfrac_l-xfrac_h)*(z-z_h)/(z_l-z_h)
		print "The xfrac cube has been interpolated using", interpolation, "interpolation."
	return xfrac

def coeval_dens(dens_dir, z, interpolation='linear'):
	"""
	dens_dir      : Give the path that contains the density-files.
	z	      : Redshift.
	interpolation : This is used when the coveal cube at that redshift is not available.
	"""
	if not interpolation in ['linear']: #, 'step', 'sigmoid', 'step_cell'
		raise ValueError('Unknown interpolation type: %s' % interpolation)
	dens_files  = glob.glob(dens_dir + '/*n_all.dat')
	dens_zs  = None
	dens_zs  = c2t.lightcone._get_file_redshifts(dens_zs, dens_files)
	if z in dens_zs:
		dens, dtype = c2t.helper_functions.get_data_and_type(dens_files[np.argwhere(z==dens_zs)])
	else:
		z_l = dens_zs[dens_zs<z].max()
		z_h = dens_zs[dens_zs>z].min()
		dens_l, dtype = c2t.helper_functions.get_data_and_type(dens_files[dens_zs[dens_zs<z].argmax()])
		dens_h, dtype = c2t.helper_functions.get_data_and_type(dens_files[dens_zs[dens_zs>z].argmin()])
		dens = dens_h + (dens_l-dens_h)*(z-z_h)/(z_l-z_h)
		print "The density cube has been interpolated using", interpolation, "interpolation."
	return dens

def coeval_vel(dens_dir, vel_dir, z, interpolation='linear'):
	"""
	vel_dir       : Give the path that contains the velocity-files.
	z	      : Redshift.
	interpolation : This is used when the coveal cube at that redshift is not available.
	"""
	if not interpolation in ['linear']: #, 'step', 'sigmoid', 'step_cell'
		raise ValueError('Unknown interpolation type: %s' % interpolation)
	dens_files = glob.glob(dens_dir + '/*n_all.dat')
	vel_files  = glob.glob(vel_dir + '/*v_all.dat')
	vel_zs     = None
	vel_zs     = c2t.lightcone._get_file_redshifts(vel_zs, vel_files)
	def get_vel(vel_file, dens_file):
		dfile = c2t.density_file.DensityFile(dens_file)
		vel_file = c2t.vel_file.VelocityFile(vel_file)
		vel = vel_file.get_kms_from_density(dfile)
		return vel
	if z in vel_zs:
		vel = get_vel(vel_files[np.argwhere(z==vel_zs)],dens_files[np.argwhere(z==vel_zs)])
	else:
		z_l = vel_zs[vel_zs<z].max()
		z_h = vel_zs[vel_zs>z].min()
		vel_l = get_vel(vel_files[vel_zs[vel_zs<z].argmax()],dens_files[vel_zs[vel_zs<z].argmax()])
		vel_h = get_vel(vel_files[vel_zs[vel_zs>z].argmin()],dens_files[vel_zs[vel_zs>z].argmin()])
		vel = vel_h + (vel_l-vel_h)*(z-z_h)/(z_l-z_h)
		print "The velocity cube has been interpolated using", interpolation, "interpolation."
	return vel

def get_adjacent_coevals(xfrac_dir, z, depth_mpc='sim'):
	"""
	xfrac_dir : Give the path that contains the xfrac-files.
	z	  : redshift.
	depth_mpc : Comoving distance to consider as dz about the provided z (Default:'sim').
		    The value given should be int/float. 'sim' takes the value from "set_sim_constants".
	"""
	xfrac_files = glob.glob(xfrac_dir + '/xfrac3d_*.bin')
	xfrac_zs = None
	xfrac_zs = c2t.lightcone._get_file_redshifts(xfrac_zs, xfrac_files)
	z_l, z_h = 0, 0
	if depth_mpc:
		if depth_mpc == 'sim': depth_mpc = c2t.conv.LB
		zl = c2t.cdist_to_z(c2t.z_to_cdist(z)-depth_mpc/2.)
		zh = c2t.cdist_to_z(c2t.z_to_cdist(z)+depth_mpc/2.)
		z_l = xfrac_zs[xfrac_zs<zl].max()
		z_h = xfrac_zs[xfrac_zs>zh].min()
	else:
		z_l = xfrac_zs[xfrac_zs<z].max()
		z_h = xfrac_zs[xfrac_zs>z].min()
	return z_l, z_h

def get_zs_list(xfrac_dir, file_type='/xfrac3d_*.bin'):
	"""
	xfrac_dir: Provide the directory whic contains all the data with redshift values in the name.
	file_type: Give the filename used to list it in Unix.
		   Example- xfrac files:   '/xfrac3d_*.bin' (Default)
			    density files: '/*n_all.dat'
	return   
	------
	The list of redshift values in the increasing order.
	"""
	xfrac_files = glob.glob(xfrac_dir + file_type)	
	xfrac_zs = None
	xfrac_zs = c2t.lightcone._get_file_redshifts(xfrac_zs, xfrac_files)
	return np.sort(xfrac_zs)

def get_dz_from_dcdist(z, dcdist=None):
	if not dcdist: dcdist = c2t.conv.LB
	zl = c2t.cdist_to_z(c2t.z_to_cdist(z)-dcdist/2.)
	zh = c2t.cdist_to_z(c2t.z_to_cdist(z)+dcdist/2.)
	return zh-zl

		 

