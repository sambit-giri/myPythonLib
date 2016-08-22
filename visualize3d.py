import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def scatter_plot(array, label, s=20, c='r'):
	"""
	Paramaters
	----------
	array (numpy array): The 3D numpy array that is to be visualized.
	label (int or numpy array of int): The labels that are to be plotted.
	s (int): The size of the plotted points (Default: 20).
	c (string): The color of the plotted points (Default: 'r').
	"""
	fig = plt.figure('3D scatter plot')
	ax = fig.add_subplot(111, projection='3d')
	if (type(label) == int):
		wh = np.argwhere(array == label)
		ax.scatter(wh[:,0], wh[:,1], wh[:,2], c=c, s=s)
		
	elif (type(label).__module__ == np.__name__):
		color = 0
		for l in label:
			wh = np.argwhere(array == l)
			ax.scatter(wh[:,0], wh[:,1], wh[:,2], c=10**color*np.ones(wh.shape[0]), s=s, cmap='jet')
			color += 1
	else: print "Enter int or numpy array as label"	
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')

	plt.show()
