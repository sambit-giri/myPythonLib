import math
import numpy as np
from random import randint
from scipy import interpolate
from interp import interp3Darray
import sys

def ray2d (ar, xth, iterations = 10000000):
	info    = ar.shape	
	lon	= max([info[0], info[1]])
	longest = math.sqrt(2*lon*lon)
	longest = math.ceil(longest)
	siz     = np.zeros(longest)
	n       = iterations
	mx      = 0                               #output gives maxmum size value

	ix       = np.arange(info[0])
	iy       = np.arange(info[1])
	ixx, iyy = np.meshgrid(ix, iy)
	iz       = ar[ixx, iyy]
	int_f    = interpolate.interp2d(ix, iy, iz, kind='cubic')


	for k in range(1, n):
		x     = randint(0, info[0]-1)
		y     = randint(0, info[1]-1)
		theta = randint(0, 360)
		
		if (ar[x,y] != 0):
			l       = math.cos(theta*3.142/180)
			m       = math.sin(theta*3.142/180)
			sz      = 0
			control = 1
			centre  = ar[x,y]
			while (control == 1):
				x = x+l
				y = y+m
				if (x < 0) or (x >= info[0]): control = 0
				if (y < 0) or (y >= info[1]): control = 0
				point = int_f(x,y)				
				if math.fabs(point - centre) > xth: control = 0 
				else: sz = sz + 1 
			
			print sz
			if sz > mx: mx = sz
			siz[sz] = siz[sz] + 1
	
	nor = np.arange(longest)
	pp  = nor*siz
	return (pp, siz, mx)


def ray3d (arr, xth, iterations = 10000000, verbose=True):
	#3D interpolation is required
	#A trilinear interpolation is used which in defined in interp.py
	
	info    = arr.shape
	lon	= max([info[0], info[1], info[2]])	
	longest = math.sqrt(3*lon*lon)
	longest = math.ceil(longest)
	siz     = np.zeros(longest)
	count   = 0

	ar  = np.zeros(arr.shape)
	ar[arr >= xth] = 1

	while(count<iterations):
		x     = randint(0, info[0]-1)
		y     = randint(0, info[1]-1)
		z     = randint(0, info[2]-1)
		if (ar[x,y,z]):
			theta = randint(0, 360)
			phi   = randint(0, 360)
			l       = math.sin(theta*3.142/180)*math.cos(phi*3.142/180)
			n       = math.sin(theta*3.142/180)*math.sin(phi*3.142/180)
			m       = math.cos(theta*3.142/180)
			sz      = 0
			control = True
			centre  = 1	
			count  += 1
			if verbose:
				perc = (count+1)*100/iterations
				msg  = str(perc) + '%'
				loading_verbose(msg)
			while (control):
				x = x+l
				y = y+m
				z = z+n
				if (x < 0) or (x >= info[0]): control = False
				if (y < 0) or (y >= info[1]): control = False
				if (z < 0) or (z >= info[2]): control = False
				point = interp3Darray(ar, x, y, z)				
				if math.fabs(point - centre) > 0.5: control = False
				else: sz += 1 
			siz[sz] = siz[sz] + 1
	
	return siz, longest

def loading_verbose(string):
	msg = ("Completed: " + string )
	sys.stdout.write('\r'+msg)
	sys.stdout.flush()


