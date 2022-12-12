import scipy.sparse.linalg as linalg
import scipy.sparse as sprs
import os
import config as conf
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
from quadtree import QuadTree, Point, Rect
from matplotlib.collections import PatchCollection

def stitch_and_solve_matrix():
    rowarray = []
    idxes = []
    colarray = []
    boundary_array = []
    M_vec = []



    data_array = [1]*len(rowarray)
    print(len(boundary_array), len(rowarray), len(colarray))
    sparse_test = sprs.coo_matrix((data_array, (rowarray, colarray)))
    csr_conversion_test = sparse_test.tocsr()
    x_test = linalg.lsqr(csr_conversion_test, M_vec)

    test_rectangles(boundary_array, x_test)

def test_rectangles(rect_array, color_vec):
    fig = plt.figure()
    ax = plt.subplot()
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    print(len(rect_array))
    patches = []
    if False:
        fig = plt.figure()
        ax = plt.subplot()
        for rect in rect_array:
            rect.draw(ax)
        plt.show()
        return
    for rect in rect_array:
        xy = (rect.west_edge, rect.north_edge)
        width = rect.east_edge - rect.west_edge
        height = rect.south_edge - rect.north_edge
        patches.append(ptc.Rectangle(xy, width, height))
    collection = PatchCollection(patches, cmap = "gray")
    collection.set_array(color_vec)
    ax.add_collection(collection)
    plt.tight_layout()
    plt.savefig("./disambiguated.png", dpi = 1440)
    print("here")
    plt.show()

def check_ranges(axis_list):
    for axis in axis_list:
        slons = []
        slats = []
        nlons = []
        nlats = []
        #print(f'axis: {axis}')
        list_of_files = os.listdir(f'{conf.point_path}a{axis}/')
        list_of_files.sort()
        for file in list_of_files:
            #print(f'file: {file}')
            with open(f'{conf.point_path}a{axis}/{file}', 'r') as fp:
                for line in fp.readlines():
                    line = line.partition("#")[0].rstrip()
                    if line:
                        npoint, spoint = line.split()
                        npoint = npoint.split(",")
                        spoint = spoint.split(",")
                        try:
                            nlons.append(float(npoint[1]))
                        except:
                            print(npoint[0])
                        try:
                            nlats.append(float(npoint[0]))
                        except:
                            print(npoint[1])
                        try:
                            slons.append(float(spoint[1]))
                        except:
                            print(spoint[0])
                        try:
                            slats.append(float(spoint[0]))
                        except:
                            print(spoint[1])
        plt.plot(nlats, nlons, '.')
        plt.plot(slats, slons, '.')
        plt.show()

def create_initial_rectfile(axis_list):
    pointlist = []
    for axis in axis_list:
        list_of_files = os.listdir(f'{conf.point_path}a{axis}/')
        list_of_files.sort()
        for file in list_of_files:
            with open(f'{conf.point_path}a{axis}/{file}', 'r') as fp:
                for line in fp.readlines():
                    line = line.partition("#")[0].rstrip()
                    if line:
                        for pair in line.split():
                            x, y = pair.split(",")
                            pointlist.append(Point(float(x), float(y)))
    print(len(pointlist))
    domain = Rect(0, 0, 360.01, 180.01)
    tree = QuadTree(domain, 16)
    for point in pointlist:
        tree.insert(point)
    boundary_list = []
    get_boundary_list(tree, boundary_list)
    print(len(boundary_list))

    path = f'{conf.matrix_path}/test_boundaries.dat'
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fp:
        fp.write("# boundary_array\n")
        for line in boundary_list:
            fp.write(f'{line}\n')


def get_boundary_list(tree, boundary_array):
    _get_boundary_list(tree, boundary_array)


def _get_boundary_list(tree, boundary_array):
    if tree.divided:
        _get_boundary_list(tree.ne, boundary_array)
        _get_boundary_list(tree.nw, boundary_array)
        _get_boundary_list(tree.sw, boundary_array)
        _get_boundary_list(tree.se, boundary_array)
    else:
        boundary_array.append(tree.boundary)


axislist = [0, 1, 2, 3, 4]

