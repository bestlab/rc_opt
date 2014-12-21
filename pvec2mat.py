#!/usr/bin/env python

import sys,math

def get_ij(ind):
	i = int(math.floor((1.0+math.sqrt(1.0+8.0*ind))/2.0))
	j = ind - i*(i-1)/2
	return i,j


for line in sys.stdin.readlines():
	ls = line.split()
	ind, wij, rij = int(ls[0]),float(ls[1]),float(ls[2])
	i,j = get_ij(ind)
	sys.stdout.write("%5i %5i %12.6e %12.6f\n" % (i,j,wij,rij))

