#!/usr/bin/env python

#convert contact list to pvec format

import sys,math

pfile = open(sys.argv[1])

def get_ij(ind):
	i = int(math.floor((1.0+math.sqrt(1.0+8.0*ind))/2.0));
	j = ind - i*(i-1)/2;
	return i+1,j+1

pdat = map(lambda x: x.split(), pfile.readlines())

for p in pdat:
	ind, w, r = int(p[0]),float(p[1]),float(p[2])
	i,j = get_ij(ind)
	print i,j,r,w

