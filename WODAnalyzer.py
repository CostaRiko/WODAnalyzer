import numpy as np
import pandas as pd
import math

class Atom:
    def __init__(self, pdb_atom_line):
        """ Read a PDB atom record and initiate Atom object:
                fields:
                    id, name, res, resSeq, x, y, z"""
        self.id = int(pdb_atom_line[7:12].replace(' ', ''))
        self.name = pdb_atom_line[13:17].replace(' ', '')
        self.res = pdb_atom_line[17:21].replace(' ', '')
        self.resSeq = int(pdb_atom_line[23:27].replace(' ', ''))
        self.x = float(pdb_atom_line[31:38].replace(' ', ''))
        self.y = float(pdb_atom_line[38: 46].replace(' ', ''))
        self.z = float(pdb_atom_line[46:55].replace(' ', ''))

    def __str__(self):
        return str(self.id)+" : "+self.name+" : "+self.res+" : "+str(self.resSeq)+\
               " : "+str(self.x)+" : "+str(self.y)+" : "+str(self.z)

class WODAnalyzer:

    rotation_sequence = []
    current_sequence = {
        'CA': None,
        'CB': None,
        'N': None
    }
    hemi_pi = math.pi/2
    protein = []
    protein_complete = False

    def __init__(self, path=None): # работаем тут
        if not path is None:
            f = open(path, 'r')

            for line in f:
                if line[0:6].find('ATOM') != -1:
                    a = Atom(line)
                    if a.resSeq is 'TIP3':
                        self.protein_complete = True
                    if not self.protein_complete:
                        if a.name == 'CA':
                            self.current_sequence['CA'] = np.array([ a.x, a.y, a.z ])
                            self.check_or_append_cs()
                        if a.name == 'CB' or a.name == 'HA1':
                            self.current_sequence['CB'] = np.array([ a.x, a.y, a.z ])
                            self.check_or_append_cs()
                        if a.name == 'N':
                            self.current_sequence['N'] = np.array([ a.x, a.y, a.z ])
                            self.check_or_append_cs()

    def check_or_append_cs(self):
        if not self.current_sequence['CA'] is None \
                and not self.current_sequence['CB'] is None \
                and not self.current_sequence['N'] is None:
            self.protein.append(self.current_sequence)
            self.clear_current_sequence()

    def clear_current_sequence(self):
        self.current_sequence = {
        'CA': None,
        'CB': None,
        'N': None
    }

    def shift(self, shift_val, vector):
        return np.array([
            vector[0] - shift_val[0],
            vector[1] - shift_val[1],
            vector[2] - shift_val[2]
        ])

    def cos_a(self, a, b):
        p1 = ( a[0] * b[0] ) + ( a[1] * b[1] ) + ( a[2] * b[2] )
        p2 = [ (a[0]**2 + a[1]**2 + a[2]**2)*(b[0]**2 + b[1]**2 + b[2]**2) ]
        return p1/np.sqrt(p2)[0]

    def x_rotation_matrix(self, sina, cosa):
        return np.matrix([
            [1, 0, 0],
            [0, cosa, -sina],
            [0, sina, cosa]
        ])

    def y_rotation_matrix(self, sina, cosa):
        return np.matrix([
            [cosa, 0, sina],
            [0, 1, 0],
            [-sina, 0, cosa]
        ])

    def z_rotation_matrix(self, sina, cosa):
        return np.matrix([
            [cosa, -sina, 0],
            [sina, cosa, 0],
            [0, 0, 1]
        ])

    def check_if_null(self, val):
        if abs(val) < 10*10**-10:
            return 0
        else:
            return val

    def clear_rotation_trigonometrical_sequence(self):
        self.rotation_sequence.clear()

    def distance(self, a, b):
        return np.sqrt([((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)])[0]

    def sin_from_cos(self, cosa):
        return np.sqrt([(1 - cosa ** 2)])[0]

    def sphere_coordinates_sequence(self, CA, CB, N):
        x_direction = np.array([1, 0, 0])
        y_direction = np.array([0, 1, 0])

        CB1 = self.shift(CA, CB)
        N1 = self.shift(CA, N)

        cb_oxy_proj = np.array([CB1[0], CB1[1], 0])
        cosa = self.check_if_null(self.cos_a(cb_oxy_proj, x_direction))
        if CB1[0] < 0:
            alpha = math.acos(cosa)
            cosa = self.check_if_null(math.cos( (self.hemi_pi - alpha) + self.hemi_pi )*-1)
        sina = self.check_if_null(self.sin_from_cos(cosa))
        if cb_oxy_proj[1] < 0:
            sina = -sina
        self.rotation_sequence.append({
            'sina': sina,
            'cosa': cosa
        })
        rmz = self.z_rotation_matrix(sina, cosa)
        CB2 = np.transpose(CB1).dot(rmz)
        N2 = np.transpose(N1).dot(rmz)
        CB2[0, 1] = self.check_if_null(CB2[0, 1])
        CB2 = np.array([CB2[0, 0], CB2[0, 1], CB2[0, 2]])
        N2 = np.array([N2[0, 0], N2[0, 1], N2[0, 2]])

        cosa = self.check_if_null(self.cos_a(CB2, x_direction))
        if CB2[0] < 0:
            alpha = math.acos(cosa)
            cosa = self.check_if_null(math.cos( (self.hemi_pi - alpha) + self.hemi_pi )*-1)
        sina = self.check_if_null(self.sin_from_cos(cosa))
        if CB2[2] > 0:
            sina = -sina
        self.rotation_sequence.append({
            'sina': sina,
            'cosa': cosa
        })
        rmy = self.y_rotation_matrix(sina, cosa)
        CB3 = np.transpose(CB2).dot(rmy)
        N3 = np.transpose(N2).dot(rmy)
        N3 = np.array([N3[0, 0], N3[0, 1], N3[0, 2]])
        CB3[0, 2] = self.check_if_null(CB3[0, 2])
        if CB3[0, 1] != 0 or CB3[0, 2] != 0:
            print("CALCULATION ERROR!")
            print(CB1, CB2, CB3)
            print(N1, N2, N3)
            print(self.rotation_sequence)
            raise IOError("Calculation error. Y and Z dimension of resulted CB vector not null. Program cancel")
        CB3 = np.array([self.check_if_null(CB3[0, 0]), self.check_if_null(CB3[0, 1]), self.check_if_null(CB3[0, 2])])
        if N3[2] != 0:
            n_yoz_proj = np.array([0, N3[1], N3[2]])
            cosa = self.cos_a(n_yoz_proj, y_direction)
            if N3[1] > 0:
                alpha = math.acos(cosa)
                cosa = self.check_if_null(math.cos((self.hemi_pi - alpha) + self.hemi_pi) * -1)
            sina = self.sin_from_cos(cosa)
            if N3[2] < 0:
                sina = -sina
            xrm = self.x_rotation_matrix(sina, cosa)
            N4 = np.transpose(N3).dot(xrm)
            N4 = np.array([self.check_if_null(N4[0, 0]), self.check_if_null(N4[0, 1]), self.check_if_null(N4[0, 2])])
            if N4[2] != 0:
                print(N1, N2, N3)
                print(self.rotation_sequence)
                raise IOError("Calculation error. Z dimension of resulted N vector not null. Program cancel")
            if N4[1] < 0:
                print(N1, N2, N3)
                print(self.rotation_sequence)
                raise IOError("Calculation error. Wrong rotational direction. N.Y cooredinate is negative. Program cancel")
        return CB3, N4


path = "C:\\Users\\softc\\Science\\poly_amino_acids\\polygly\\frames\\frames.pdb"
w = WODAnalyzer(path)
CA = np.array([1, 1, 1])
CB = np.array([-2, -2, -3])
N = np.array([-3, -2, -2])
print(w.sphere_coordinates_sequence(CA, CB, N))
