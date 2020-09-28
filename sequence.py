from pdb import *

class Sequence:

    def __init__(self, pdb=""):
        if pdb != "":
             self.pdb = PDB(pdb).parse()

