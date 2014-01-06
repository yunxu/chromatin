#!/usr/bin/python
import sys
import os
import warnings
import Bio.PDB,numpy

from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.PDBIO import PDBIO
from numpy import *

#--------------------------------------------------------------------
if len(sys.argv) != 3:
	print >>sys.stderr, "Usage: {0} ref.pts out.pdf".format(sys.argv[0])
	sys.exit()

ref_ptsfilename, outputfilename = sys.argv[1:4]

if not os.path.exists(ref_ptsfilename) :
	sys.stderr.write("ERROR: pts file %r was not found!\n" % (ref_ptsfilename))
	sys.exit()

# --------------------------------------------------------------------
FNameIndex = "../script/data/analysis/ENm008/ENm008_GM12878_Cont_Index.txt"
IndexListOld = [line.strip() for line in open(FNameIndex)]
# index should start from 0
IndexList = [int(x)-1 for x in IndexListOld]

scale=0.01

#--------------------------------------------------------------------
def ReadXYZ	(filename,scale=1):
	""" Read uniformly distributed sphere points from file
	"""
	num_lines = sum(1 for line in open(filename))
	points = zeros(shape=(num_lines,3))
	ind = 0
	for line in open(filename).readlines():
		points[ind] = array((line.split()[0:3])) 
		points[ind] = points[ind] * scale
		ind = ind+1
	return points	

#--------------------------------------------------------------------
# ref_ptsfilename = "K562.pts"
refid = "ref"
structure = Structure(refid)
model_ref = Model(1)
chain_ref = Chain("A")
points_ref = ReadXYZ(ref_ptsfilename,scale)
	
num_count = 0
for i in range(0,shape(points_ref[IndexList])[0]):
	num_count = num_count +1
	# res_id = (' ',IndexList[i],' ')
	res_id = (' ',num_count,' ')	
	residue = Residue(res_id,'ALA',' ')
	cur_coord = tuple(points_ref[IndexList[i]])
	# atom = Atom('CA',cur_coord,0,0,' ',num_count,num_count,'C')
	atom = Atom('CA',cur_coord,0,0,' ',num_count,num_count,'C')
	residue.add(atom)
	chain_ref.add(residue)
model_ref.add(chain_ref)
structure.add(model_ref)

#--------------------------------------------------------------------
print ("scale=%s" % scale)

#--------------------------------------------------------------------
io=PDBIO()
io.set_structure(structure)
fn = outputfilename
io.save(fn)
fout = open(fn,"a")
for i in range(1,shape(points_ref[IndexList])[0]):
	fout.write( "CONECT%5d%5d\n" % (i, i+1))
print "output file: " + outputfilename

#--------------------------------------------------------------------
for i in range(0,len(IndexList)-1):
	dist = numpy.linalg.norm(points_ref[IndexList[i]] - points_ref[IndexList[i+1]])
	print ("%d %d %.3f" % (IndexList[i]+1,IndexList[i+1]+1,dist))

#--------------------------------------------------------------------

