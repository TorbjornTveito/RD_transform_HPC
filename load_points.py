import os

import config as conf
from quadtree import Rect, Point



'''
Sugested way of loading data from files

Data is saved under ./data/[a1 - ax]/[000001 - 999999].dat
format is X1:float,Y1:float X2:float,Y2:float M:float # Comments and blank lines are ok.
Note the space ' ' between point1 and point2

This function will process each file/line/point individually, limiting the amount of memory needed

'''
def get_all_points():
    idx = 0
    for axis in conf.matrix_axes:
        list_of_files = os.listdir(f'{conf.path}a{axis}/')
        list_of_files.sort()
        for file in list_of_files:
            with open(f'{conf.path}a{axis}/{file}', 'r') as fp:
                for line in fp.readlines():
                    line = line.partition("#")[0].rstrip()
                    if line:
                        idx += 1
                        for p in line.split():
                            p = tuple(map(float, p.split(',')))
                            yield Point(*p, idx)



'''
This function will only return points that are within the given rect
'''
def get_points(rect):
    for point in get_all_points():
        if rect.contains(point):
            yield point



'''
Read all rects from file conf.rect_file
stored in file as:
cx:float,cy:float,w:float,h:float # Comments and blank lines are ok
'''
def get_all_rects():
    with open(conf.rect_file, 'r') as fp:
        rects = (row.partition("#")[0].rstrip() for row in fp)
        rects = (row for row in rects if row)
        rects = (tuple(map(float, row.split(','))) for row in rects)
        rects = [Rect(*rect_tuple) for rect_tuple in rects]

    return rects


'''
Possibly not needed...
'''
def get_rect(num, rects):
    return rects[num]



def main():

    points = get_points(get_rect(0))
    [points]
    print(list(points))

if __name__ == "__main__":
    main()