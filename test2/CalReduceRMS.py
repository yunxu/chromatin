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
if len(sys.argv) != 4:
	print >>sys.stderr, "Usage: {0} ref.pts alt.pts out.pdf".format(sys.argv[0])
	sys.exit()

ref_ptsfilename, alt_ptsfilename, outputfilename = sys.argv[1:4]

if not os.path.exists(ref_ptsfilename) or not os.path.exists(alt_ptsfilename):
	sys.stderr.write("ERROR: pts file %r or %r was not found!\n" % (ref_ptsfilename, alt_ptsfilename))
	sys.exit()

# --------------------------------------------------------------------

FNameIndex = "ENm008_GM12878_Cont_Index.txt"
IndexListOld = [line.strip() for line in open(FNameIndex)]
# index should start from 0
IndexList = [int(x)-1 for x in IndexListOld]

scale=.01

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
ref_ptsfilename = "K562.pts"
refid = "ref"
structure = Structure(refid)
model_ref = Model(1)
chain_ref = Chain("A")
points_ref = ReadXYZ(ref_ptsfilename,scale)
	
num_count = 0
for i in range(0,shape(points_ref[IndexList])[0]):
	num_count = num_count +1
	res_id = (' ',num_count,' ')
	residue = Residue(res_id,'ALA',' ')
	cur_coord = tuple(points_ref[IndexList[i]])
	atom = Atom('CA',cur_coord,0,0,' ',num_count,num_count,'C')
	residue.add(atom)
	chain_ref.add(residue)
model_ref.add(chain_ref)
structure.add(model_ref)

#--------------------------------------------------------------------
altid = "alt"
structure_alt = Structure(refid)
model_alt = Model(2)
chain_alt = Chain("A")
points_alt = ReadXYZ(alt_ptsfilename,scale)
	
num_count = 0
for i in range(0,shape(points_alt[IndexList])[0]):
	num_count = num_count +1
	res_id = (' ',num_count,' ')
	residue = Residue(res_id,'ALA',' ')
	cur_coord = tuple(points_alt[IndexList[i]])
	atom = Atom('CA',cur_coord,0,0,' ',num_count,num_count,'C')
	residue.add(atom)
	chain_alt.add(residue)
model_alt.add(chain_alt)
structure.add(model_alt)

#--------------------------------------------------------------------
super_imposer = Bio.PDB.Superimposer()
ref_atoms = list(model_ref.get_atoms())
alt_atoms = list(model_alt.get_atoms())
super_imposer.set_atoms(ref_atoms, alt_atoms)
super_imposer.apply(model_alt.get_atoms())
print ("%s %s RMS= %.3f" % (ref_ptsfilename,alt_ptsfilename,super_imposer.rms))

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


