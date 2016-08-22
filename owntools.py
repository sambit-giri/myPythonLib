import numpy as np
import c2raytools as c2t
from scipy.signal import argrelextrema
from sklearn.cluster import KMeans

def threshold_21cm(cube, p):
	return (1.-p)*cube.min() + p*cube.max()

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

def threshold_kmeans(cube, upper_lim=False):
	"""
	The input is the brightness temperature cube.
	"""
	array = np.zeros(cube.shape)
	km = KMeans(n_clusters=2)
	X  = cube.reshape(-1,1)
	y = KMeans(n_clusters=2).fit_predict(X)
	if X[y==1].mean()>X[y==0].mean(): t_th = X[y==0].max()
	else: t_th = X[y==1].max()
	if upper_lim: array[cube<=t_th] = 1
	else: array[cube>t_th] = 1
	print "The output contains a tuple with binary-cube and determined-threshold."
	return array, t_th
	

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
	dummy = zz[zz!=0].min()/10.
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
					if i*i+j*j+k*k <= radius*radius: array[x,y,z] = label
	if periodic: print "Periodic sphere of radius", radius, "made at", centre
	else: print "Non-periodic sphere of radius", radius, "made at", centre
				
				
					 

