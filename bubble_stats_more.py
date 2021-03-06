from sklearn.cluster import DBSCAN
import numpy as np
import datetime, time
from Friends_of_Friends import FoF_search
from scipy import ndimage as ndi

from skimage.morphology import watershed
from skimage.feature import peak_local_max

def fof(data, xth=0.5):
	"""
	FOF bubble
	
	Parameter
	---------
	input  : 3D array of ionization fraction.
	xth    : The threshold value (Default: 0.5).

	Output
	------
	The output is a tuple containing output-map and volume-list array.
	"""
	t1 = datetime.datetime.now()
	out_map, size_list = FoF_search(data, xth)
	t2 = datetime.datetime.now()
	runtime = (t2-t1).total_seconds()/60

	print "Program runtime: %f minutes." %runtime
	print "The output is a tuple containing output-map and volume-list array respectively."

	return out_map, size_list


def dbscan_cube(cube, eps=1.12, min_samples=7, metric='euclidean', weight=True, upper_lim=True):
	"""
	DBSCAN for bubble size
	
	Parameter
	---------
	input      : 3D array of ionization fraction, i.e., range of elements is [0,1].
	eps        : Minimum distance from neighbours (Default:1.12).
	min_samples: Minimum number of neighbours including the point itself (Default: 7).
	metric     : Define the distance metric (Default: 'euclidean').
		     - from scikit-learn: ['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan']
                     - from scipy.spatial.distance: ['braycurtis', 'canberra', 'chebyshev',
                                          'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski',
                                          'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto',
                                          'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath',
                                          'sqeuclidean', 'yule']
	weight     : It is a boolean defining if each point should be weighted with the point value
		     as the cluster forming probability (Default: False).

	Output
	------
	The output is a tuple containing output-map and volume-list array.
	"""
	t1 = datetime.datetime.now()
	if upper_lim: cube = -1.*cube
	mn, mx = cube.min(), cube.max()
	X = np.zeros((cube.shape[0]*cube.shape[1]*cube.shape[2],4))
	X[:,:3] = np.argwhere(~np.isnan(cube))
	for i in xrange(cube.shape[0]*cube.shape[1]*cube.shape[2]):
		X[i,3] = cube[X[i,0],X[i,1],X[i,2]]
		if not weight: X[i,3] = (X[i,3]-mn)*cube.shape[0]/(mx-mn)
		else: X[i,3] = (X[i,3]-mn)*min_samples/(mx-mn)
	if weight: db = DBSCAN(eps=eps, min_samples=min_samples, metric=metric).fit(X[:,:3],sample_weight=X[:,-1])
	else: db = DBSCAN(eps=eps, min_samples=min_samples, metric=metric).fit(X)
	labels = db.labels_
	print "The number of clusters formed is %d." %labels.max()
	if labels.max()<=0:
		print "Try other values of eps and min_samples."
		return 0, 0

	out_cube = np.zeros(cube.shape)
	sizes = np.zeros(labels.max())
	for i in xrange(cube.shape[0]*cube.shape[1]*cube.shape[0]):
		out_cube[X[i,0],X[i,1],X[i,2]] = labels[i]+1
	for j in xrange(len(sizes)):
		sizes[j] = labels[labels==j+1].shape[0]

	t2 = datetime.datetime.now()
	runtime = (t2-t1).total_seconds()/60
	print "Program runtime: %f minutes." %runtime
	return out_cube, sizes

def watershed_cube(data, xth=0.5):
	"""
	Watershed
	
	Parameter
	---------
	input  : 3D array of ionization fraction.
	xth    : The threshold value (Default: 0.5).

	Output
	------
	The output is a tuple containing output-map and volume-list array.
	"""
	image = data>xth
	distance = ndi.distance_transform_edt(image)
	local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((3, 3, 3)), labels=image)
	markers = ndi.label(local_maxi)[0]
	labels = watershed(-distance, markers, mask=image)
	sizes  = [labels[labels==i].size for i in np.unique(labels)]
	return labels, sizes

def dist_from_volumes(sizes, resolution=1., bins=100):
	"""
	Volume distribution and effective radius distribution.
	
	Parameter
	---------
	sizes      : List of volumes in pixel**3 units.
	resolution : Distance between two pixels in cMpc (Default: 1).
	bins       : Number of bins of volumes in log space.

	Output
	------
	The output is a tuple conatining 4 numpy arrays: V, VdP/dV, r, rdp/dr.
	The radius calculated here is the effective radius.
	"""
	vols  = np.array(sizes)
	radii = (vols*3./4./np.pi)**(1./3.)
	ht_v  = np.histogram(np.log10(vols), bins=bins)
	ht_r  = np.histogram(np.log10(radii), bins=bins)
	vs, d_v  = np.zeros(len(ht_v[0])+1), np.zeros(len(ht_v[0])+1)
	vs       = 10.**ht_v[1]*resolution**3
	d_v[:-1] = 1.*ht_v[0]/np.sum(ht_v[0])
	dummy = d_v[d_v!=0].min()/1000.
	d_v[d_v==0] = dummy
	rs, d_r  = np.zeros(len(ht_r[0])+1), np.zeros(len(ht_r[0])+1)
	rs       = 10.**ht_r[1]*resolution
	d_r[1:]  = 1.*ht_r[0]/np.sum(ht_r[0])
	d_r[0]   = d_r[d_r!=0].min()/1000.
	print "The output is a tuple conatining 4 numpy array: V, VdP/dV, r, rdp/dr."
	return vs, d_v, rs, d_r
	

def get_distribution(array, resolution=1., bins=100, sizes=False):
	if sizes:
		sizes = array
	else:
		mn, mx = array.min(), array.max()
		sizes  = np.arange(mx)+1.
		for i in xrange(int(mx)):
			label = i+1
			sizes[i] = len(array[array==label])
			print label,
	ht   = np.histogram(np.log(sizes), bins=bins)
	vols, dist = np.zeros(len(ht[0])+1), np.zeros(len(ht[0])+1)
	vols      = np.exp(ht[1])*resolution
	dist[:-1] = ht[0]

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




