#!/usr/bin/python
import sys
import warnings
import Bio.PDB,numpy
import argparse

from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.PDBIO import PDBIO
from numpy import *

#--------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Read pts file, out pdb file, according to scale')
parser.add_argument('--src', required=True,  help='pts file')
parser.add_argument('--dst', required=False, help='pdb file')
parser.add_argument('--scale', required=False, default='1.0', type=float, help='default 1')
parser.add_argument('--link', required=False, default=False, action='store_true', help='default false')
parser.add_argument('--bfactor', required=False, help='bfactor file')
parser.add_argument('--column', required=False, default=1,type=int, help='bfactor column selection')

args = vars(parser.parse_args())
#--------------------------------------------------------------------
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
#--------------------------------------------------------------------
def ReadBfactor(filename,column):
    """ Read bfactor file
    """
    num_lines = sum(1 for line in open(filename) if not line.startswith('#'))
    bfactors = zeros(shape=(num_lines,1))
    ind = 0
    for line in open(filename).readlines():
        if not line.startswith('#'):
            bfactors[ind] = array((line.split())[column]) 
            ind = ind+1
    return bfactors	
#--------------------------------------------------------------------
points = ReadXYZ ( args['src'], args['scale'])
if ( args['bfactor'] is not None):
    print "read bfactor file column %d" % args['column']
    bfactors = ReadBfactor(args['bfactor'],args['column'])
else:
    bfactors = zeros(len(points))

model = Model(1)
chain = Chain("A")
structure = Structure("ref")

num_count = 0
for i in range(0,shape(points)[0]):
    num_count = num_count +1
    res_id = (' ',num_count,' ')
    residue = Residue(res_id,'ALA',' ')
    cur_coord = tuple(points[i])
    bfactor = bfactors[i]
    atom = Atom('CA',cur_coord,bfactor,0,' ','CA',num_count,'C')
    residue.add(atom)
    chain.add(residue)

model.add(chain)
structure.add(model)
# --------------------------------------------------------------------
io=PDBIO()
io.set_structure(structure)
if ( args['dst'] is None):
    fn = sys.stdout
    io.save(fn)
    if ( args['link'] ):
        for i in range(1,shape(points)[0]):
            fn.write( "CONECT%5d%5d\n" % (i, i+1))
else:
    fn = args['dst']
    io.save(fn)
    fout = open(fn,"a")
    if (args['link'] ):
        for i in range(1,shape(points)[0]):
            fout.write( "CONECT%5d%5d\n" % (i, i+1))
    fout.close()
#    print "output file: " + args['dst']


