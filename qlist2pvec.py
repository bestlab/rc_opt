#!/usr/bin/env python

#convert contact list to pvec format

import sys,math

qfile = open(sys.argv[1])

qlist = []
wsum = 0.0
for q in qfile.readlines():
	qs = q.split()
	i,j,rij = int(qs[0])-1,int(qs[1])-1,float(qs[2])
	idx = j*(j-1)/2+i
	if len(qs) == 4:
		wij = float(qs[3])
		wsum += abs(wij)
		qlist.append((idx,rij,wij))
	else:
		qlist.append((idx,rij))

ww = 1.0/float(len(qlist))

for q in qlist:
	if len(q) == 2:
		wij = ww
		idx,rij=q
	else:
		idx,rij,w=q
		wij = w/wsum
	sys.stdout.write("%5i %12.6e %12.6e\n"%(idx,wij,rij))




