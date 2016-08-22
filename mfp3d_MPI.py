import numpy as np
import sys
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE
import c2raytools as c2t
from math import *
from random import randint
from interp import interp3Darray


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


cube = np.load(sys.argv[1]+'.npy')
flname = str(sys.argv[2])
xth = float(sys.argv[3])
n    = int(sys.argv[4])
verbose = str(sys.argv[5])

ar = np.zeros(cube.shape)
ar[cube>=xth] = 1
info = ar.shape

local_n = n/size
count = 0
siz = []

while(count < local_n):
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
		if rank == size: 
				if verbose == 'True': print "Nodes used", size, "| iter", count*size
				elif verbose != 'False': print verbose, "| Nodes used", size, "| iter", count*size
		siz.append(sz)
	
if rank == 0:
        total = np.array(siz)
        for i in range(1, size):
		recv_buffer = np.empty(local_n, dtype=np.int)
                comm.Recv(recv_buffer, ANY_SOURCE)
                total = np.hstack((total, recv_buffer))
else:
        # all other process send their result
        comm.Send(np.array(siz))

print "Process ", rank, "is done"

if comm.rank == 0:
	np.savetxt(flname, total)

	
