from math import *
import numpy as np
from random import randint
from scipy import interpolate
import sys

def mfp2d(arr, xth=0.5, iterations=1000000, verbose=True):
	info    = arr.shape	
	lon	= max([info[0], info[1]])
	longest = sqrt(2*lon*lon)
	longest = ceil(longest)
	siz     = np.zeros(longest)
	count   = 0

	ar  = np.zeros(arr.shape)
	ar[arr >= xth] = 1

	while(count<iterations):
		x     = randint(0, info[0]-1)
		y     = randint(0, info[1]-1)
		if verbose:
			perc = (k+1)*100/iterations
			msg  = str(perc) + '%'
			loading_verbose(msg)
		
		if (ar[x,y]):
			theta = randint(0, 360)
			l       = cos(theta*3.142/180)
			m       = sin(theta*3.142/180)
			sz      = 0
			control = True
			centre  = 1
			count  += 1
			while (control == 1):
				x = x+l
				y = y+m
				if (x < 0) or (x >= info[0]): control = False
				if (y < 0) or (y >= info[1]): control = False
				point = interp2Darray(ar, x, y, z)				
				if fabs(point - centre) > 0.5: control = False
				else: sz = sz + 1 
			siz[sz] = siz[sz] + 1
	size_px = np.arange(longest)
	return siz, size_px


def mfp3d(arr, xth=0.5, iterations=10000000, verbose=True):
	#3D interpolation is required
	#A trilinear interpolation is used which in defined in interp.py
	
	info    = arr.shape
	lon	= max([info[0], info[1], info[2]])	
	longest = sqrt(3*lon*lon)
	longest = ceil(longest)
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
			l       = sin(theta*3.142/180)*cos(phi*3.142/180)
			n       = sin(theta*3.142/180)*sin(phi*3.142/180)
			m       = cos(theta*3.142/180)
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
				if (z < 0) or (z >= info[1]): control = False
				point = interp3Darray(ar, x, y, z)				
				if fabs(point - centre) > 0.5: control = False
				else: sz += 1 
			siz[sz] = siz[sz] + 1
	size_px = np.arange(longest)
	return siz, size_px

def interp2Darray(sq, x, y, missing = 0):
	#Biliniear interpolation for periodic lattice
	#Missing point considered zero
	missing = 0

	info = sq.shape
	l = info[0]
	m = info[1]

	x0 = floor(x)
	x1 = x0+1
	y0 = floor(y)
	y1 = y0+1	

	xd = (x-x0)/(x1-x0)
	yd = (y-y0)/(y1-y0)

	V00 = 'empty'
	V01 = 'empty'
	V10 = 'empty'
	V11 = 'empty'

	if (x0 < 0 or x0 >= l):
		V00 = missing
		V01 = missing
		
	if (x1 < 0 or x1 >= l):
		V10 = missing
		V11 = missing

	if (y0 < 0 or y0 >= m):
		V00 = missing
		V10 = missing

	if (y1 < 0 or y1 >= m):
		V01 = missing
		V11 = missing


	if V00 == 'empty': V00 = sq[x0,y0]
	if V01 == 'empty': V01 = sq[x0,y0]
	if V10 == 'empty': V10 = sq[x0,y1]
	if V11 == 'empty': V11 = sq[x0,y1]


	c0 = V00*(1-xd) + V10*xd	
	c1 = V01*(1-xd) + V11*xd

	c = c0*(1-yd)+c1*yd

	return c



def interp3Darray(cube, x, y, z, missing = 0):
	#Triliniear interpolation for periodic lattice
	#Missing point considered zero
	missing = 0

	info = cube.shape
	l = info[0]
	m = info[1]
	n = info[2]

	x0 = floor(x)
	x1 = x0+1
	y0 = floor(y)
	y1 = y0+1	
	z0 = floor(z)
	z1 = z0+1

	xd = (x-x0)/(x1-x0)
	yd = (y-y0)/(y1-y0)
	zd = (z-z0)/(z1-z0)

	V000 = 'empty'
	V001 = 'empty'
	V010 = 'empty'
	V011 = 'empty'
	V100 = 'empty'
	V101 = 'empty'
	V110 = 'empty'
	V111 = 'empty'

	if (x0 < 0 or x0 >= l):
		V000 = missing
		V001 = missing
		V010 = missing
		V011 = missing

	if (x1 < 0 or x1 >= l):
		V100 = missing
		V101 = missing
		V110 = missing
		V111 = missing

	if (y0 < 0 or y0 >= m):
		V000 = missing
		V001 = missing
		V100 = missing
		V101 = missing

	if (y1 < 0 or y1 >= m):
		V010 = missing
		V011 = missing
		V110 = missing
		V111 = missing

	if (z0 < 0 or z0 >= n):
		V000 = missing
		V010 = missing
		V100 = missing
		V110 = missing

	if (z1 < 0 or z1 >= n):
		V001 = missing
		V011 = missing
		V101 = missing
		V111 = missing

	if V000 == 'empty': V000 = cube[x0,y0,z0]
	if V001 == 'empty': V001 = cube[x0,y0,z1]
	if V010 == 'empty': V010 = cube[x0,y1,z0]
	if V011 == 'empty': V011 = cube[x0,y1,z1]
	if V100 == 'empty': V100 = cube[x1,y0,z0]
	if V101 == 'empty': V101 = cube[x1,y0,z1]
	if V110 == 'empty': V110 = cube[x1,y1,z0]
	if V111 == 'empty': V111 = cube[x1,y1,z1]

	c00 = V000*(1-xd) + V100*xd
	c10 = V010*(1-xd) + V110*xd
	c01 = V001*(1-xd) + V101*xd
	c11 = V011*(1-xd) + V111*xd

	c0 = c00*(1-yd) + c10*yd	
	c1 = c01*(1-yd) + c11*yd

	c = c0*(1-zd)+c1*zd

	return c

def loading_verbose(string):
	msg = ("Completed: " + string )
	sys.stdout.write('\r'+msg)
	sys.stdout.flush()

