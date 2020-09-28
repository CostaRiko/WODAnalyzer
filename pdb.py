
class Atom:
    def __init__(self, pdb_atom_line):
        """ Read a PDB atom record and initiate Atom object:
                fields:
                    id, name, res, resSeq, x, y, z"""
        self.id = int(pdb_atom_line[7:12].replace(' ', ''))
        self.name = pdb_atom_line[13:17].replace(' ', '')
        self.res = pdb_atom_line[17:21].replace(' ', '')
        self.resSeq = int(pdb_atom_line[23:27].replace(' ', ''))
        self.coordinates = {
            'x' : float(pdb_atom_line[31:38].replace(' ', '')),
            'y' : float(pdb_atom_line[38:46].replace(' ', '')),
            'z' : float(pdb_atom_line[46:55].replace(' ', ''))
        }

    def __str__(self):
        return str(self.id)+" : "+self.name+" : "+self.res+" : "+str(self.resSeq)+\
               " : "+str(self.x)+" : "+str(self.y)+" : "+str(self.z)

class Residue:
    def __init__(self, s):
        self.seq = s
        self.CA = s.get_by_atomname('CA')
        self.CB = s.get_by_atomname('CB')
        self.N = s.get_by_atomname('CN')

class Sequence(object):
    pos = 0
    hoh = []
    def __init__(self, atom_sequence):
        self.atoms = atom_sequence

    def size(self):
        return len(self.atoms)

    def __getitem__(self, item):
        return self.atoms[item]

    def __iter__(self):
        return self

    def __next__(self):
        self.pos += 1
        if self.pos >= len(self.atoms):
            self.pos = 0
            raise StopIteration
        else:
            return self.atoms[self.pos]

    def __str__(self):
        line = ""
        for x in range(len(self.atoms)):
            line += "\t:\t"+self.atoms[x].name
            if x%30 == 0 and x != 0:
                line += '\n'
        return line

    def get_by_atomname(self, a_name_condition):
        CAs = []
        for a in self.atoms:
            if a.name == a_name_condition:
                CAs.append(a)
        return Sequence(CAs)

    def get_by_res_id(self, res_id):
        res = []
        for a in self.atoms:
            if a.resSeq == res_id:
                res.append(a)
        return Sequence(res)

class PDB:
    atoms = []
    hoh = []
    sequences = []
    residues = []
    max_res_id = 0
    def __init__(self, path=""):
        if path == "":
            print("Error of PDB input")
            exit()
        f = open(path, 'r')
        self.content = f.readlines()
        f.close()

    def show_pdb_text_line(self):
        for l in self.content:
            print(l)
        print("sizeof(l) = {}".format(len(self.content)))

    def parse(self):
        for record in self.content:
            if record[0: 4] == 'ATOM':
                a = Atom(record)
                if a.resSeq > self.max_res_id:
                    self.max_res_id = a.resSeq
                self.atoms.append(a)
            if record[0: 3] == 'TER':
                self.sequences.append(Sequence(self.atoms))
                self.atoms = []
        print("chains: {}".format(len(self.sequences)))

    def parse_first_model(self):
        first_model_complete = False
        for record in self.content:
            if record[0: 4] == 'ATOM':
                a = Atom(record)
                if a.res != "TIP3":
                    if a.resSeq > self.max_res_id:
                        self.max_res_id = a.resSeq
                    self.atoms.append(a)
            if record[0:3] == 'END':
                self.sequences.append(Sequence(self.atoms))
                self.atoms = []
                break

    def parse_water(self):
        current_id = -1
        hoh_mol = []
        for record in self.content:
            if record[0: 4] == 'ATOM':
                a = Atom(record)
                if a.res == "TIP3":
                    if a.resSeq != current_id:
                        if current_id == -1:
                            current_id = a.resSeq
                            self.hoh.append(hoh_mol)
                            continue
                        current_id = a.resSeq
                        self.hoh.append(hoh_mol)
                        hoh_mol = []
                        hoh_mol.append(a)
                    else:
                        hoh_mol.append(a)

class CalculatedObject:
    CA = None
    CB = None
    N = None
    HOH = []

    def __init__(self, ca, cb, n, hoh):
        self.CA = ca
        self.CB = cb
        self.N = n
        self.HOH = hoh

