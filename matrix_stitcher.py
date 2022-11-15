import scipy.sparse.linalg as linalg
import scipy.sparse as sprs
import config as conf
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
from matplotlib.collections import PatchCollection

def stitch_and_solve_matrix():
    rowarray = []
    idxes = []
    colarray = []
    boundary_array = []
    filecount = len(os.listdir(conf.output_dir))
    for i in range(filecount):
        rowarray.append(load_pickle_file("/matrix/rowarray_" + str(i) + ".pkl"))
        idxes.append(load_pickle_file("/matrix/idxes" + str(i) + ".pkl"))
        colarray.append(load_pickle_file("/matrix/colarray" + str(i) + ".pkl"))
        boundary_array.append(load_pickle_file("/matrix/boundary_array" + str(i) + ".pkl"))

    print(len(boundary_array), len(rowarray_test), len(colarray_test))
    sparse_test = sprs.coo_matrix((data_array, (rowarray_test, colarray_test)))
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