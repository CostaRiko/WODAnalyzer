from pdb import *
import numpy as np

def distance(a, b):
    return np.sqrt([( (a.x - b.x)**2 + (a.y - b.y)**2 + (a.z - b.z)**2 ) ])[0]

def check_hydrogen_bounds(a, b):
    if a.name.find('H') != -1 and a.name.find('N') != -1 and a.name.find('O') != -1:
        return True
    else:
        return False

def encount_H_bounds(aa, bb, i, ii):
    H_bounds_count = 0
    chains_count = []
    for a in aa.atoms:
        for b in bb.atoms:
            d = distance(a, b)
            if d < 4:
                H_bounds_count += 1
    print("H-bonds between chain {} and chain {} is {}".format(i, ii, H_bounds_count))



test_file = "C:\\Users\\softc\\Science\\2q8b.pdb"
p = PDB(test_file)
p.parse()

seq = p.sequences
for i, a in enumerate(seq):
    seq.remove(a)
    for ii, b in enumerate(seq):
        encount_H_bounds(a, b, i, ii)