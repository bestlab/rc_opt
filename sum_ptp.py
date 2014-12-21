#!/usr/bin/env python

import sys

def read_ptp(file):
	data = {}
	for line in open(file).readlines():
		if line[0] == "#":
			ls = line.split()
			if len(ls) == 10:
				ptp = float(ls[9])
			elif len(ls) > 2 and ls[2] == 'length':
				w = float(ls[7])
			elif len(ls) > 1 and ls[1] == 'Number':
				ntp = int(ls[6])
			continue
		q, pq, pqtp, ptpq = map(float, line.split())
		data[q] = [ pq, pqtp, ptpq ]
	return data, ptp, w, ntp

data = []
ptp = []
weights = []
NTP = 0
nf = len(sys.argv) - 1
for f in sys.argv[1:]:
	x = read_ptp(f)
	data.append(x[0])
	ptp.append(x[1])
	weights.append(x[2])
	NTP += x[3]

#weights = [ 4.0, 4.0, 1.0, 1.0, 1.0, 1.0 ]

qs = data[0].keys()
qs.sort()

pTP = 0.0
ptpqdat = {}
wsum = 0.0
for i in range(nf):
	#w = weights[i]
	w = 1.
	wsum += w
	pTP += ptp[i]*w
	for q in qs:
		pq, pqtp = data[i][q][0], data[i][q][1]
		if q not in ptpqdat.keys():
			ptpqdat[q] = [ 0.0, 0.0 ]
		ptpqdat[q][0] += pq*w;
		ptpqdat[q][1] += pqtp*w;

pTP /= wsum
sys.stdout.write("# pTP = %12.5f\n" % (pTP))
sys.stdout.write("# NTP = %i\n" % (NTP))

for q in qs:
	if ptpqdat[q][0] != 0.0:
		pTPq = ptpqdat[q][1]*pTP/ptpqdat[q][0]
	else:
		pTPq = 0.0
	sys.stdout.write("%12.5f %12.5f %12.5f %12.5f\n" % (q,
			ptpqdat[q][0]/wsum, ptpqdat[q][1]/wsum, pTPq))
