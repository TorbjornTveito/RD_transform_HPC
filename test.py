import matplotlib.pyplot as plt
from array import array
from quadtree import Point
import struct
import config as conf

import os



def load_points():
    with open('./data/point_data/a0/3_points.dat', 'r') as fp:
        data = [line for line in fp]
    p1 = []
    p2 = []
    for line in data:
        if line[0]== '#':
            continue

        points = line.split(' ')
        p1x, p1y = points[0].split(',')
        p2x, p2y = points[1].split(',')
        p1.append((float(p1x), float(p1y)))
        p2.append((float(p2x), float(p2y)))

    return p1, p2

def save_bin_points():
    #
    p1 = Point(43.5256, 123.5948, 1) ###1
    p2 = Point(47.1324, 149.3433, 2) ###2
    p3 = Point(51.4832, 134.4383, 3) ###3
    p4 = Point(-45.4596, 68.5324, 1) ###1
    p5 = Point(-38.4921, 28.4382, 2) ###2
    p6 = Point(-28.2483, 48.8384, 3) ###3

    bin = struct.pack('ffff', p1.x, p1.y, p4.x, p4.y)

    p2.x, p2.y, p5.x, p5.y = struct.unpack('ffff', bin)

    print(p1)
    print(p2)

    print(p4)
    print(p5)



def load_bin_points(rect):
    with open('./data/point_data/a0/8_points.bin', 'rb') as fp:
        data = array('d')
        data.fromstring(fp.read())
        print(data)
        data =data.tolist()
        print(data)

        data = list(zip(data[0::2], data[1::2]))

        for i, point in enumerate(data):
            if rect.contains(point):
                yield (Point(*point, i//2))






def plot_points():
    points1, points2 = load_points()

    points1_lat = [point[0] for point in points1]
    points1_lon = [point[1] for point in points1]
    points2_lat = [point[0] for point in points2]
    points2_lon = [point[1] for point in points2]

    plt.plot(points1_lon, points1_lat, '.',  c = "red")
    plt.plot(points2_lon, points2_lat, '.',  c = "blue")
    plt.show()

#load_bin_points()

def convert(axis_list):
    idx = 0
    for axis in axis_list:
        #print(f'axis: {axis}')
        list_of_files = os.listdir(f'{conf.point_path}a{axis}/')
        list_of_files.sort()
        for file in list_of_files:
            if file[-4:] == '.dat':
                new_file = file[:-4] + '.bin'
                pointlist = []

                with open(f'{conf.point_path}a{axis}/{file}', 'r') as fp:
                    for line in fp.readlines():
                        line = line.partition("#")[0].rstrip()
                        if line:
                            idx += 1
                            for p in line.split():
                                pointlist.extend(p.split(','))

                pointlist = [float(p) for p in pointlist]
                pointlist = array('d', pointlist)
                with open(f'{conf.point_path}a{axis}/{new_file}', 'wb') as fp:
                    pointlist.tofile(fp)
axislist = [0, 1, 2, 3, 4]
convert(axislist)





