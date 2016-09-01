from sklearn.cluster import DBSCAN
import numpy as np
import datetime, time
from Friends_of_Friends import FOF_search

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
	out_map, size_list = FOF_search(data, xth)
	t2 = datetime.datetime.now()
	runtime = (t2-t1).total_seconds()/60

	print "Program runtime: %f minutes." %runtime
	print "The output is a tuple containing output-map and volume-list array respectively."

	return out_map, size_list

def dbscan_cube(cube, eps=1.12, min_samples=7, metric='auto', weight=False):
	"""
	DBSCAN for bubble size
	
	Parameter
	---------
	input      : 3D array of ionization fraction.
	eps        : Minimum distance from neighbours (Default:1.12).
	min_samples: Minimum number of neighbours including the point itself (Default: 7).
	metric     : Define the distance metric (Default: 'auto').
	weight     : It is a boolean defining if each point should be weighted with the point value
		     as the cluster forming probability (Default: False).

	Output
	------
	The output is a tuple containing output-map and volume-list array.
	"""
	mn, mx = cube.min(), cube.max()
	X = np.zeros((cube.shape[0]*cube.shape[1]*cube.shape[2],4))
	X[:,:3] = np.argwhere(~np.isnan(cube))
	for i in xrange(cube.shape[0]*cube.shape[1]*cube.shape[2]):
		X[i,3] = cube[X[i,0],X[i,1],X[i,2]]
		X[i,3] = (X[i,3]-mn)*cube.shape[0]/(mx-mn)
	if weight: db = DBSCAN(eps=eps, min_samples=min_samples, metric=metric).fit(X,sample_weight=X[:,-1]*(mx-mn)/cube.shape[0])
	else: db = DBSCAN(eps=eps, min_samples=min_samples, metric=metric).fit(X)
	labels = db.labels_
	out_cube = np.zeros(cube.shape)
	for i in xrange(cube.shape[0]*cube.shape[1]*cube.shape[0]):
		out_cube[X[i,0],X[i,1],X[i,2]] = labels[i]+1
	for j in xrange(len(sizes)):
		sizes[j] = labels[labels==j+1].shape[0]
	print "Program runtime: %f minutes." %runtime
	return out_cube, sizes
