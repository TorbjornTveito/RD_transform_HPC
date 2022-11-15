from quadtree import Point, Rect, QuadTree
import numpy as np
import os

import config as conf

from mpi4py import MPI

RUN_WITHOUT_MPI = True

comm = MPI.COMM_WORLD
mpi_rank = comm.Get_rank()
mpi_size = comm.Get_size()
mpi_admin = mpi_size - 1


'''
Get the current rectangle to populate based on the rectangle index
'''
def get_rect(num):
    with open(conf.rect_file, 'r') as fp:
        rects = (row.partition("#")[0].rstrip() for row in fp)
        rects = (row for row in rects if row)
        rects = [tuple(map(float, row.split(','))) for row in rects]
        rect = Rect(*rects[num])
    return rect

def get_num_rects():
    with (conf.rect_file, 'r') as fp:
        rects = (row.partition("#")[0].rstrip() for row in fp)
        rects = [row for row in rects if row]

    return len(rects)

'''
Get all points contained inside rect in all axes that are specified in the axis_list
'''
def get_points(axis_list, rect):
    idx = 0
    for axis in axis_list:
        #print(f'axis: {axis}')
        list_of_files = os.listdir(f'{conf.point_path}a{axis}/')
        list_of_files.sort()
        for file in list_of_files:
            #print(f'file: {file}')
            with open(f'{conf.point_path}a{axis}/{file}', 'r') as fp:
                for line in fp.readlines():
                    line = line.partition("#")[0].rstrip()
                    if line:
                        idx += 1
                        for p in line.split():
                            p = tuple(map(float, p.split(',')))
                            if rect.contains(p):
                                yield Point(*p, idx)


'''
Do the magic
'''
def calculate(num):
    rect = get_rect(num)
    rect = Rect(0, 0, 360.01, 180.01)
    grid_tree = QuadTree(rect, conf.point_limit)

    colarray = []
    rowarray = []
    for point in get_points(conf.pixel_axes, rect):
        grid_tree.insert(point)
    boundary_array = []
    get_boundary_list(grid_tree, boundary_array)
    num_boundaries = len(boundary_array)

    for i, rect in enumerate(boundary_array):
        print(num_boundaries-i)
        idxes = [point.payload for point in get_points(conf.matrix_axes, rect)]
        rowarray += idxes
        colarray += [i]*len(idxes)


    write_output(rowarray, colarray, boundary_array, num)


def write_output(rowarray, colarray, boundary_array, job_num):
    path = f'{conf.matrix_path}/{job_num}_rows.dat'
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fp:
        fp.write("# Rowarray\n")
        for line in rowarray:
            fp.write(f'{line}\n')

    path = f'{conf.matrix_path}/{job_num}_cols.dat'
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fp:
        fp.write("# colarray\n")
        for line in colarray:
            fp.write(f'{line}\n')

    path = f'{conf.matrix_path}/{job_num}_boundaries.dat'
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fp:
        fp.write("# boundary_array\n")
        for line in boundary_array:
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




def run():
    # Start by calculating the job corresponding to mpi_rank
    num = mpi_rank
    while True:
        # Calculate current job
        calculate(num)

        # Signal admin that a new job can be deligated to this worker
        comm.send(mpi_rank, dest=mpi_admin, tag=1)
        num = comm.recv(mpi_admin, tag=2)

        # Shut down if all jobs are complete
        if num == -1:
            break


def run_admin():
    # Distribute workloads to workers until all jobs has completed
    for i in range(mpi_admin, get_num_rects()):
        worker = comm.recv(source=MPI.ANY_SOURCE, tag=1)
        comm.send(i, dest=worker, tag=2)

    # Signal all workers that the work is complete
    for i in range(0, mpi_size -1):
        worker = comm.recv(source=MPI.ANY_SOURCE, tag=1)
        comm.send(-1, dest=worker, tag=2)

def main():
    if (mpi_rank == mpi_admin):
        # MPI admin will distribute jobs to workers
        run_admin()
    else:
        # Workers will cary out work
        run(mpi_rank)

def main_local_test():
    num = 0
    calculate(num)


if __name__ == "__main__":
    if RUN_WITHOUT_MPI:
        main_local_test()
    else:
        main()
