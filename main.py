import pandas as pd, os
from PIL import Image
import numpy as np
import sys

path = "C:\\Users\\softc\\Science\\poly_amino_acids\\"+sys.argv[1]+"\\frames\\frames.pdb"
log = pd.read_csv(path+'.distribution.log', delimiter="\t", names=['zenit', 'azimut', 'ZO', 'AO', 'RO'], dtype=float)

dimensions = (3840, 2160)

def convert_radians_to_pixels(rad, range='x'):
    global dimensions
    if range == 'x':
        point = round( float((rad/(3.14))*dimensions[1]))-1
    elif range == 'y':
        point = round( float((rad/(3.14)) * dimensions[0]))-1
    return point

def convert_radius_to_color(radius):
    return 100+int(radius*155/50+((3.14**5)*2))


print(log)
count = 0
ITERATIONS = 10000000
lenlog = len(log)
step = int(lenlog/ITERATIONS)
print(log['zenit'].max(), log['azimut'].max(), log['zenit'].min(), log['azimut'].min())
print('POLYGLY')
print('LENLOG: {}'.format(lenlog))
a = np.zeros((dimensions[1], dimensions[0], 3), dtype=np.uint8)
for i, aa in enumerate(log.iloc):
    x = convert_radians_to_pixels(float(aa['zenit']))
    y = convert_radians_to_pixels(float(aa['azimut']), range='y')
    #print(x, aa.zenit)
    #print(y, aa.azimut)
    if aa['RO'] > 0 and aa['RO'] <= 5:
        a[x][y][0] = 255
        a[x][y][1] = 0
        a[x][y][2] = 0
    if aa['RO'] > 5 and aa['RO'] <= 10:
        a[x][y][0] = 0
        a[x][y][1] = 255
        a[x][y][2] = 0
    if aa['RO'] > 10 and aa['RO'] <= 15:
        a[x][y][0] = 0
        a[x][y][1] = 0
        a[x][y][2] = 255
    if a[x][y][1] != 0:
        a[x][y][1] += 1
    else:
        a[x][y][1] = 100
    print("PROGRAM STATE: index={}, completion_percents={}".format(i, i*100/lenlog), end='\r')
    #if i % ITERATIONS == 0 and i != 0:
    #    break
#for i in range(100):
#    y = i
#    for x in range(255):
#        val = 100 + int(x * 155 / 500 + ((3.14 ** 5) * 2))
#        a[y][x][0] = 0
#        a[y][x][1] = 255 - x
#        a[y][x][2] = 255 - x

print('saving figure')
i = Image.fromarray(a)
i.save(path+'.distribution.log.test'+".png")
i.show()
#plt.savefig(path+'.distribution.log'+".png")
