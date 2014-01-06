#!/usr/bin/python
import sys
import os
import warnings

from numpy import *

# number_sample_points = 64

#--------------------------------------------------------------------
if len(sys.argv) != 3:
	print >>sys.stderr, "Usage: {0} test.pts test.dist".format(sys.argv[0])
	sys.exit()

ptsfilename, outputfilename = sys.argv[1:3]
collisionlength = 100
collisionlength = 300
if not os.path.exists(ptsfilename):
	sys.stderr.write("ERROR: pts file %r was not found!\n" % (ptsfilename))
	sys.exit()
# --------------------------------------------------------------------
scale=.01
scale=1

# --------------------------------------------------------------------
def ReadXYZ	(filename,scale=1):
	""" Read uniformly distributed sphere points from file
	"""
	num_lines = sum(1 for line in open(filename) if not line.startswith('#'))
	points = zeros(shape=(num_lines,3))
	ind = 0
	for line in open(filename).readlines():
		if not line.startswith('#'):
			points[ind] = array((line.split())[0:3]) 
			points[ind] = points[ind] * scale
			ind = ind+1
	return points	
	
# --------------------------------------------------------------------
def inner_prod(v1, v2):
     'inner production of two vectors.'
     sum = 0
     for i in xrange(len(v1)):
            sum += v1[i] * v2[i]
     return sum

# --------------------------------------------------------------------


num_nodes = 1
points = ReadXYZ(ptsfilename,scale)

num_count = 0
num_row = shape(points)[0]
num_col = shape(points)[1]

distance_tuples = []
for i in range(0,num_row-1):
    for j in range(i+1,num_row):
        delta_point = points[i] - points[j]
        distance = sqrt(inner_prod(delta_point, delta_point))
        if distance < collisionlength:
            distance_tuples.append((i,j,distance))

distance_tuples_sorted = sorted(distance_tuples,key=lambda distance: distance[2])

f = open(outputfilename,'w')
for i in range(0,len(distance_tuples_sorted)):
    row = distance_tuples_sorted[i][0]+1
    col = distance_tuples_sorted[i][1]+1
    distance = distance_tuples_sorted[i][2]
    print >> f, repr(row).rjust(6), repr(col).rjust(6), ("%0.2f" % (distance)).rjust(12) 
f.close()

# --------------------------------------------------------------------

# 
