import numpy as np
from scipy.stats import spearmanr
import c2raytools as c2t
import matplotlib.pyplot as plt
import owntools
from skimage.segmentation import slic, mark_boundaries
from skimage.filters import threshold_otsu

def slic_cube(cube, n_segments=1000, compactness=10.0, max_iter=20, sigma=0, min_size_factor=0.5, max_size_factor=3, cmap=None):
	if cmap is not None: 
		color   = plt.get_cmap(cmap)
		multichannel = True
		cube = color(cube)
		cube = np.delete(cube, 3, -1)
	else:
		multichannel = False
	labels = slic(cube, n_segments=n_segments, compactness=compactness, max_iter=max_iter, sigma=sigma, max_size_factor=max_size_factor, slic_zero=True, multichannel=multichannel)
	return labels

def see_label(out_map, label):
	binary = np.zeros(out_map.shape)
	if out_map[out_map==label].size: binary[out_map==label] = 1
	else: print "The entered label in not present in the map."
	return binary

def binary_stitch(data, labels, stitch='mean', thres=None):
	X1 = data.reshape(-1,1)
	X  = np.argwhere(data!=np.nan)*X1
	y  = labels.reshape(-1,1)
	y1 = [X1[y==i].mean() for i in np.unique(y)]
	if not thres:
		if stitch=='otsu': thres = threshold_otsu(np.array(y1))
		else: thres = X1.mean()
	y2 = np.zeros(y.shape)
	for i in np.unique(y): y2[y==i] = y1[i]
	y2 = y2 < thres
	return y2.reshape(data.shape)

def mark_boundaries(sl, lab, mark=True):
	assert sl.ndim == 2 and lab.ndim == 2
	bd  = mark_boundaries(np.zeros(sl.shape), lab)
	out = sl.copy()
	if mark: out[bd[:,:,0]==1] = sl.max()*1.01
	else: out[0,0] = sl.max()*1.01
	return out

	
