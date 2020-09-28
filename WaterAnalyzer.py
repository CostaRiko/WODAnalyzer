import threading as tr
from multiprocessing import cpu_count
import numpy as np
import pandas as pd
import time, os
from pdb import *
import math

# Выполняем многопоточный анализатор ориентационного распределения воды

# Технические задачи
# Выполнить чтение PDB

mutex = False
DISTANCE = 15
CAs_objects = []
calculated_objects_sequence = []
CB = None
N = None
calculated_object = None
hohs = []

class WaterAnalyzer:

    rotation_sequence = []

    x_direction = np.array([1, 0, 0])
    y_direction = np.array([0, 1, 0])

    csv_str = ""



    def distance(self, a, b):
        return np.sqrt([((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)])[0]

    def analyze(self, path):
        """
        Метод рассчета значений ориентационного распределения воды
        :return None:
        """
        global DISTANCE, CAs_objects, calculated_objects_sequence, CB, N, calculated_object, hohs
        p = PDB(path)
        p.parse_first_model()
        p.parse_water()
        atoms = p.sequences[0]
        hoh = p.hoh
        CAs = atoms.get_by_atomname('CA')
        if path.find('polygly') != -1:
            CBs = atoms.get_by_atomname('HA2')
        else:
            CBs = atoms.get_by_atomname('CB')
        Ns = atoms.get_by_atomname('N')

        for i, ca in enumerate(CAs):
            CB = CBs.get_by_res_id(ca.resSeq)[0]
            N = Ns.get_by_res_id(ca.resSeq)[0]
            for h in hoh:
                if self.distance(ca, h[1]) < DISTANCE:
                    hohs.append(h)
            calculated_object = CalculatedObject(ca, CB, N, hohs)
            calculated_objects_sequence.append(calculated_object)
        for i, calc_obj in enumerate(calculated_objects_sequence):
            while tr.active_count() >= cpu_count():
                print(i, tr.active_count(), cpu_count())
                time.sleep(1)
            t = tr.Thread(target= self.intermediate_calculations, args=(calc_obj, ))
            t.start()
        while tr.active_count() != 1:
            print(i, tr.active_count(), cpu_count())
            time.sleep(1)
        self.write_orientational_distribution_file(path)


    def intermediate_calculations(self, calc_obj):
        calc_obj.CB.x = calc_obj.CB.x - calc_obj.CA.x
        calc_obj.CB.y = calc_obj.CB.y - calc_obj.CA.y
        calc_obj.CB.z = calc_obj.CB.z - calc_obj.CA.z
        calc_obj.N.x = calc_obj.N.x - calc_obj.CA.x
        calc_obj.N.y = calc_obj.N.y - calc_obj.CA.y
        calc_obj.N.z = calc_obj.N.z - calc_obj.CA.z
        a = np.array([calc_obj.CB.x, calc_obj.CB.y, calc_obj.CB.z])
        b = np.array([calc_obj.N.x, calc_obj.N.y, calc_obj.N.z])
        a, b, rotation_sequence = self.run_coordinate_update_sequence(a, b)
        if a[0, 1] != 0 or a[0, 2] != 0 or b[0, 2] != 0:
            print("Вектор CB не приведен к требуемому положению. Последовательность поворотов не верна!")
            quit()
        self.update_water_orientational_stat(calc_obj.HOH, calc_obj.CA, rotation_sequence)
        print('ПОТОК ВЫЧИСЛЕНИЯ ЗАВЕРШЕН')

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
        if val < 10*10**-10:
            return 0
        else:
            return val

    def clear_rotation_trigonometrical_sequence(self):
        self.rotation_sequence.clear()

    def run_coordinate_update_sequence(self, a, b): # Векторы должны выходить из начала координат
        """
        Function is creating rotation sin and cos sequences
        :param a: CB vector
        :param b: N vector
        :return: updated_cb, updated_n
        """
        rotation_sequence = []
        a_XOY_proj = np.array([a[0], a[1], 0])
        cosa = self.cos_a(self.x_direction, a_XOY_proj)
        sina = np.sqrt( [ (1 - cosa**2) ] )[0]
        if a_XOY_proj[1] < 0:
            sina = -sina
        rotation_sequence.append({'sin': sina, 'cos': cosa})

        rmz = self.z_rotation_matrix(sina, cosa)
        a1 = np.transpose(a).dot(rmz)
        b1 = np.transpose(b).dot(rmz)
        a1[0, 1] = self.check_if_null(a1[0, 1])

        a1 = np.array([ a1[0, 0], a1[0, 1], a1[0, 2] ])
        b1 = np.array([ b1[0, 0], b1[0, 1], b1[0, 2] ])
        cosa = self.cos_a(self.x_direction, a1)
        sina = np.sqrt([(1 - cosa ** 2)])[0]
        hemi_pi = math.pi/2
        if a1[0] < 0:
            alpha = math.acos(cosa)
            cosa = math.cos((hemi_pi - alpha) + hemi_pi) * -1
        if a1[2] > 0:
            sina = -sina
        rotation_sequence.append({'sin': sina, 'cos': cosa})

        rmy = self.y_rotation_matrix(sina, cosa)
        a2 = np.transpose(a1).dot(rmy)
        b2 = np.transpose(b1).dot(rmy)
        a2[0, 2] = self.check_if_null(a2[0, 2])
        b2 = np.array([b2[0, 0], b2[0, 1], b2[0, 2]])
        b2_XOY_proj = np.array([0, b2[1], b2[2]])
        cosa = self.cos_a(b2_XOY_proj, self.y_direction)
        sina = np.sqrt([(1 - cosa ** 2)])[0]
        if a1[0] < 0:
            alpha = math.acos(cosa)
            cosa = math.cos((hemi_pi - alpha) + hemi_pi) * -1
        if b2[2] < 0:
            sina = -sina
        rotation_sequence.append({'sin': sina, 'cos': cosa})

        rmx = self.x_rotation_matrix(sina, cosa)
        b3 = np.transpose(b2).dot(rmx)
        b3[0, 2] = self.check_if_null(b3[0, 2])
        return a2, b3, rotation_sequence

    def water_coordinates_update(self, hoh, CA, rs):
        """
        Выполняет обновление координат воды согласно последовательности поворотов
        Возврашает обновленные координаты
        :param hoh:
        :return update_hoh:
        """
        for i, a in enumerate(hoh):
            vector = np.array([hoh[i].x - CA.x, hoh[i].y - CA.y, hoh[i].z  - CA.z])
            rmz = self.z_rotation_matrix(rs[0]['sin'], rs[0]['cos'])
            rmy = self.y_rotation_matrix(rs[1]['sin'], rs[1]['cos'])
            rmx = self.x_rotation_matrix(rs[2]['sin'], rs[2]['cos'])
            vector = np.transpose(vector).dot(rmz)
            vector = vector.dot(rmy)
            vector = vector.dot(rmx)
            hoh[i].x = vector[0, 0]
            hoh[i].y = vector[0, 1]
            hoh[i].z = vector[0, 2]
        return hoh

    def write_orientational_data(self, hoh):
        global mutex
        try:
            h = hoh.get_by_atomname('OH2')[0]
        except:
            return
        OH2 = hoh.get_by_atomname('OH2')[0]

        H1 = hoh.get_by_atomname('H1')[0]
        H2 = hoh.get_by_atomname('H2')[0]

        H1.x = H1.x - OH2.x
        H1.x = H1.y - OH2.y
        H1.x = H1.z - OH2.z

        H2.x = H2.x - OH2.x
        H2.y = H2.y - OH2.y
        H2.z = H2.z - OH2.z

        zenit_H1 = np.arccos( H1.z / np.sqrt( H1.x**2 + H1.y**2 + H1.z**2 ) )
        zenit_H2 = np.arccos(H2.z / np.sqrt(H2.x ** 2 + H2.y ** 2 + H2.z ** 2))

        azimut_H1 = np.arctan(H1.y / H1.x)
        azimut_H2 = np.arctan(H2.y / H2.x)

        zenit_OH2 = np.arccos( OH2.z / np.sqrt( OH2.x**2 + OH2.y**2 + OH2.z**2 ) )
        azimut_OH2 = np.arctan( OH2.y / OH2.x )
        radius_OH2 = np.sqrt( OH2.x**2 + OH2.y**2 + OH2.z**2 )

        while mutex:
            time.sleep(1)
        mutex = True
        self.csv_str += str(zenit_H1) + "\t" + str(azimut_H1) + "\t" + \
                        str(zenit_OH2) + "\t" + str(azimut_OH2) + "\t" + \
                        str(radius_OH2) + "\n"
        self.csv_str += str(zenit_H2) + "\t" + str(azimut_H2) + "\t" + \
                        str(zenit_OH2) + "\t" + str(azimut_OH2) + "\t" + \
                        str(radius_OH2) + "\n"
        mutex = False

    def update_water_orientational_stat(self, hohs, CA, rs):
        for i, hoh in enumerate(hohs):
            hohs[i] = self.water_coordinates_update(hoh, CA, rs)
            print("FIRST STEP: {}".format(i * 100 / len(hohs)), end='\r')
        for i, hoh in enumerate(hohs):
            self.write_orientational_data(Sequence(hoh))
            print("SECOND STEP: {}".format(i * 100 / len(hohs)), end='\r')

    def write_orientational_distribution_file(self, path):
        f = open(path+'.distribution.log', 'w+')
        f.write(self.csv_str)
        f.close()


path = "C:\\Users\\softc\\Science\\poly_amino_acids\\polyala\\frames\\frames.pdb"
w = WaterAnalyzer()
print('WaterAnalyzer object is initialized')
w.analyze(path)