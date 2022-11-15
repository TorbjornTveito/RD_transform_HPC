import numpy as np
import config as conf
import generator_parameters as gen
import os
import matplotlib.pyplot as plt

import PIL
import scipy.interpolate as intp


from mpi4py import MPI
from try_spicey import find_sub, find_dist
import util
import datetime as datetime

comm = MPI.COMM_WORLD
mpi_rank = comm.Get_rank()
mpi_size = comm.Get_size()
mpi_admin = mpi_size - 1

#mpi_rank = mpi_size = mpi_admin = 0


# Useful constants
R_M = 1737.4e3


# run function
def find_range(axis_num, job_num):
    #SKRIV ET PROGRAM SOM LES UT EN LISTE FRA EN FIL OG VELGE ELEMENT JOB_NUM
    pass

def find_max_doppler(rang, bw):
    max_dop = np.sqrt(2*rang/R_M - (rang/R_M)**2 ) * bw / 2
    return(max_dop)

def dist(srp_lat,srp_lon,nlat=300,nlon=300):
    latgrid=np.linspace(-90.0,90.0,num=nlat)
    longrid=np.linspace(-180.0,180.0,num=nlon)
    lats,lons=np.meshgrid(latgrid,longrid)
    x,y,z=R_M*util.lunar2cartesian(lats,lons)
    sc=util.lunar2cartesian(srp_lat,srp_lon)
    d=R_M-(x*sc[0]+y*sc[1]+z*sc[2])
    inv=np.where(d > R_M)
    d[inv]=np.nan
    d = np.ma.masked_invalid(d)
    return(d,latgrid,longrid)

def find_dop(sublat1, sublon1, sublat2, sublon2, dt, axis_params):
    d1,_,_ = dist(sublat1,sublon1)
    d2,_,_ = dist(sublat2,sublon2)

    return 2*axis_params.obs_frequency*(d2-d1)/dt.seconds/gen.c

'''
Finds the angle between the north pole and the apparent rotation axis
given Doppler geometry in the form of two sub-radar points
'''
def find_np_angle(sublat1, sublon1, sublat2, sublon2):
    coord1 = util.lunar2cartesian(sublat1, sublon1)
    coord2 = util.lunar2cartesian(sublat2, sublon2)

    xaxis = [1, 0, 0]
    yaxis = [0, 1, 0]
    zaxis = [0, 0, 1]

    ymat, zmat = aquire_matrices(coord1)

    rotated_vec_1 = np.dot(zmat, np.dot(ymat, coord1))
    rotated_vec_2 = np.dot(zmat, np.dot(ymat, coord2))

    rotated_diff = rotated_vec_1 - rotated_vec_2

    np_vec = np.cross(rotated_vec_1, rotated_diff)
    norm_np_vec = np_vec / np.sqrt(np_vec[0]**2 + np_vec[1]**2 + np_vec[2]**2)

    np_angle = np.arccos(np.dot(norm_np_vec, zaxis)) * 180 / np.pi

    crossprod = np.cross(norm_np_vec, zaxis)
    if (np.dot(xaxis, crossprod) < 0):
        np_angle = -np_angle
    return(np_angle)


'''
Finds rotation matrix to rotate a given angle around y-axis.
'''
def yrot_mat(theta):
    yrot = np.array([[np.cos(theta), 0 , np.sin(theta)], [0, 1 , 0], [-np.sin(theta), 0, np.cos(theta)]])
    return(yrot)

'''
Finds rotation matrix to rotate a given angle around z-axis.
'''
def zrot_mat(theta):
    zrot = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    return(zrot)

'''
Finds rotation matrix to rotate a given angle around x-axis.
'''
def xrot_mat(theta):
    xrot = np.array([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
    return(xrot)

'''
Finds rotation matrices .
'''
def aquire_matrices(coord_vec):
    z_angle = np.arccos(np.dot([coord_vec[0], coord_vec[1]], [1, 0])/np.linalg.norm([coord_vec[0], coord_vec[1]]))
    zmat = zrot_mat(-(z_angle))

    partial_rot = np.dot(zmat, coord_vec)

    y_angle = np.arccos((np.dot([partial_rot[0], partial_rot[2]], [1,0])/(np.linalg.norm([partial_rot[0], partial_rot[2]]))))
    ymat = yrot_mat(y_angle)

    return(ymat, zmat)

def rotate_to_srp(coord_vec, sublat, sublon): # makes SRP be 1, 0, 0 in cartesian
    if sublon > 180:
        sublon = sublon - 360

    coord_vec = np.array(coord_vec)
    yrot = yrot_mat(np.float128(sublat* np.pi/180))
    partial_rotation = np.dot(yrot, coord_vec)
    zrot = zrot_mat(np.float128(-sublon*np.pi/180))

    full_rot = np.dot(zrot, partial_rotation)

    return(full_rot)

def rotate_from_srp(coord_vec, sublat, sublon): # inverse of rotate to srp
    if sublon > 180:
        sublon = sublon - 360

    coord_vec = np.array(coord_vec)

    yrot = yrot_mat(sublat * np.pi / 180)

    partial_rotation = np.dot(coord_vec, yrot)

    zrot = zrot_mat(-sublon * np.pi / 180)

    full_rot = np.dot(partial_rotation, zrot)
    return(full_rot)

def pixel_area(lat):
    return(gen.c * gen.wavelength * gen.baud_length / (4 * gen.obsdur * gen.ANGULAR_ROT * np.abs(np.sin(np.pi*lat/180))))

def find_aoi(rang):
    x = (R_M - rang)/R_M
    y = np.sqrt(np.abs(1-x**2))
    point = [1, 0, 0]
    aoi = np.arccos(np.dot(point, [x, y, 0]))
    return(aoi)

def scattering_law(aoi, const = 1000):
    sigma = const*(1/((np.cos(aoi))**4 + const *(np.sin(aoi))**2))**(3/2)
    return(sigma)


def slaw_over_data(axis_name):
    filecount = len(os.listdir(DATA_DIR + axis_name + '/sums'))
    for i in range(filecount):
        range_vec = util.load_pickle(axis_name + "/ranges/" + str(i) + "_pickled")
        sum_vec = util.load_pickle(axis_name + "/sums/" + str(i) + "_pickled")
        if len(range_vec) == 0:
            continue
        if isinstance(range_vec, list):
            rang = range_vec[0]
        else:
            rang = range_vec
        sum_vec = remove_slaw_effects(rang, sum_vec)
        picklefunc(axis_name + "/sums/" + str(i) + "_pickled", sum_vec)

def remove_slaw_effects(rang, point_sum_list, lunar_dist):
    slaw = scattering_law(find_aoi(rang))
    sprcs = slaw * gen.gain**2 * gen.tx_power / (4*np.pi * (lunar_dist + rang)**2)**2
    #print(sprcs)
    point_sum_list = point_sum_list/sprcs
    return(point_sum_list)

def create_dops(rang, axis_params, srp1, srp2, bw, fres, mapfunc, axis_num, job_num):
    lunar_dist = find_dist(str(datetime.datetime(*axis_params.obsdate)), axis_params.radar_lat, axis_params.radar_lon, axis_params.radar_el)
    if axis_params.do_SPRCS:
        if axis_params.do_SLAW:
            aoi = find_aoi(rang)
            slaw = scattering_law(aoi)
        else:
            slaw = 1
        sprcs = slaw * axis_params.gain**2 * axis_params.tx_power / (4*np.pi * (lunar_dist + rang)**2)**2
    else:
        sprcs = 1

    max_dop = find_max_doppler(rang, bw)

    pointlist = []
    point_sums = []
    dop = 0
    k = 0
    print("number of dops: ", len(np.concatenate((np.arange(0, max_dop, fres), np.arange(-fres, -max_dop, -fres)))))
    for dop in np.concatenate((np.arange(0, max_dop, fres), np.arange(-fres, -max_dop, -fres))):
        points, point_sum = find_average_return(rang, dop, srp1, srp2, bw, sprcs, axis_params, fres, mapfunc)#(rang, dop, srp1, srp2, bw, rcs)
        if axis_params.make_SAR:
            if point_sum is None:
                continue
            point_sums.append(point_sum)
        pointlist.append(points)
    if axis_params.make_SAR:
        write_output(pointlist, axis_num, job_num, point_sums)
    else:
        write_output(pointlist, axis_num, job_num, [])


def create_point(rang, dop, srp1, srp2, bw):
    x = (R_M - rang)/R_M
    y = (-dop / (bw/2))
    z = np.sqrt(np.abs(1 - x**2 - y**2))                            #  the abs() should not actually matter,
    np_angle = find_np_angle(srp1[0], srp1[1], srp2[0], srp2[1])    #  but i think that due to binning, the length
    point = np.array([x, y, z])                                     #  of the xy vector can go over 1. which is bad.
    pair_point = np.array([x, y, -z])                               #  but its not much bigger than 1. which is good.
                                                                    #  less than 0.001 diff, i think.
    np_correction_matrix = xrot_mat(-np_angle * np.pi / 180)
    point_rotation = np.dot(np_correction_matrix, point)
    pair_rotation = np.dot(np_correction_matrix, pair_point)

    point_sub = rotate_from_srp(point_rotation, srp1[0], srp1[1])
    pair_sub = rotate_from_srp(pair_rotation, srp1[0], srp1[1])

    lat1, lon1 = util.cartesian2lunar(point_sub)
    lat2, lon2 = util.cartesian2lunar(pair_sub)
    if lon1 > 180:
        lon1 = lon1 - 360
    if lon2 > 180:
        lon2 = lon2 - 360
    effective_latitude = util.cartesian2lunar(point)[0]
    return([lon1, lat1], [lon2, lat2], effective_latitude)

def supersample(r, dop, delr, deldop,n):
    coord_list = np.linspace(1/(2*n)-0.5, 0.5-1/(2*n), n)

    r_vec = [r + delr * k for k in coord_list]
    dop_vec = [dop + deldop * k for k in coord_list]
    return(r_vec, dop_vec)



def find_average_return(rang, dop, srp1, srp2, bw, sprcs, axis_params, fres, mapfunc):
    point1, point2, elat = create_point(rang, dop, srp1, srp2, bw)
    if not axis_params.make_SAR:
        return((point1, point2), 0)

    if axis_params.do_area:
        if elat > 0.01:
            area = pixel_area(elat)
        else:
            return((0, 0), None)
    else:
        area = 1


    samples = 5
    super_r, super_dop = supersample(rang, dop, axis_params.range_res, fres, samples)
    north_real = []
    north_im = []
    south_real = []
    south_im = []
    for ridx in range(samples):
        for didx in range(samples):
            point1, point2, elat = create_point(super_r[ridx], super_dop[didx], srp1, srp2, bw)

            north = (mapfunc(point1[1], point1[0])+0.01) # surface scattering simulation
            south = (mapfunc(point2[1], point2[0])+0.01)

            north_adjusted = north # spherical backscatter + ray propagation
            south_adjusted = south # multiply by pixel area for total power

            north_real.append(np.random.normal(0, np.sqrt(north_adjusted / 2))) # creating real voltage value
            north_im.append(1j*np.random.normal(0, np.sqrt(north_adjusted / 2))) # creating imaginary voltage value

            south_real.append(np.random.normal(0, np.sqrt(south_adjusted / 2))) # these are gaussians with variance equal to signal power
            south_im.append(1j*np.random.normal(0, np.sqrt(south_adjusted / 2)))

    north_real = np.array(north_real)
    north_im = np.array(north_im)
    south_real = np.array(south_real)
    south_im = np.array(south_im)

    total_real = north_real + south_real
    total_im = north_im + south_im

    total_power = area * sprcs * ((total_real + total_im) * (total_real - total_im)).real

    if False:
        print("area: ", area)
        print("sprcs: ", sprcs)
        print("north", north)

    point_sum = np.mean(total_power).real

    return((point1, point2), point_sum)

def find_srps(axisparams, dt = 60):
    obsdate1 = datetime.datetime(*axisparams.obsdate)
    srp1 = find_sub(str(obsdate1), axisparams.radar_lat, axisparams.radar_lon, axisparams.radar_el)
    obsdate2 = obsdate1 + datetime.timedelta(seconds = dt)
    srp2 = find_sub(str(obsdate2), axisparams.radar_lat, axisparams.radar_lon, axisparams.radar_el)
    return(srp1, srp2)


def load_map():
    lats_interp = np.linspace(-90, 90, 512)
    lons_interp = np.linspace(-180, 180, 1024)
    im = np.flip(np.array(PIL.Image.open('./data/generator_input/lroc3.jpg').convert('L')).T, axis = 1)
    mapfunc = intp.interp2d(lats_interp, lons_interp, im)
    return(mapfunc)




'''
må ha:
range fil, srp1 fil og srp2 fil, LROC bilde, axis folder, range resolution, frequency resolution.
'''



def write_output(points, axis_num, job_num, point_sums):
    path = f'{conf.point_path}/a{axis_num}/{job_num}_points.dat'
    os.makedirs(os.path.dirname(path), exist_ok=True)
    data = [f'{p[0][0]},{p[0][1]} {p[1][0]},{p[1][1]}' for p in points]
    if len(points) == 0:
        return()
    with open(path, 'w') as fp:
        fp.write("# P1 P2\n")
        for line in data:
            fp.write(f'{line}\n')

    path = f'{conf.measurement_path}/a{axis_num}/{job_num}_meas.dat'
    with open(path, 'w') as fp:
        fp.write("# M\n")
        for line in point_sums:
            fp.write(f'{line}\n')


def create_job_list():
    job_list = []
    for axis_num, axis in enumerate(gen.list_of_axes):
        srp1, srp2 = find_srps(axis)
        num_jobs = R_M / axis.range_res
        for job_num in range(num_jobs):
            job_list.append((axis_num, job_num, srp1, srp2))
    return job_list


def run_admin():
    job_list = create_job_list()
    for job in job_list:
        w = comm.recv(source=MPI.ANY_SOURCE, tag=1)
        comm.send(job, dest=w, tag=2)
    for w in range(mpi_size - 1):
        w = comm.recv(source=MPI.ANY_SOURCE, tag=1)
        comm.send(-1, dest=w, tag=2)


def execute_job(job, mapfunc):
    axis_num = job[0]
    job_num = job[1]
    srp1 = job[2]
    srp2 = job[3]
    axis_params = gen.list_of_axes[axis_num]
    dops = find_dop(srp1[1], srp1[0], srp2[1], srp2[0], datetime.timedelta(seconds = 60), axis_params)
    bw = np.max(dops)*4 * np.pi
    fres = np.pi * 2 * axis_params.RD_decimation_factor * 1 / axis_params.obsdur
    if axis_params.range_res * job_num <= R_M:
        rang = axis_params.range_res*job_num
    else:
        print("fucked range. fix?!")
        quit()
    create_dops(rang, axis_params, srp1, srp2, bw, fres, mapfunc, axis_num, job_num)


def run():
    mapfunc = load_map()
    while True:
        comm.send(mpi_rank, dest=mpi_admin, tag=1)
        job = comm.recv(source=mpi_admin, tag=2)
        if job == -1:
            break

        execute_job(job, mapfunc)

def run_test():
#    print_sub_to_check(0)
    print(len(gen.list_of_axes))
    for axis_num, axis in enumerate(gen.decimated_list_of_axes):
        print(f'Calculating axis {axis_num}')
        axis_params = axis
        srp1, srp2 = find_srps(axis_params)
        print(srp1, srp2)
        mapfunc = load_map()
        jobs = []
        rang = 0
        k = 0
        while rang <= R_M:
            jobs.append(k)
            k +=1
            rang = k*axis_params.range_res
        for i in jobs:
            execute_job((axis_num, i, srp1, srp2), mapfunc)


def main():
    if (mpi_rank == mpi_admin):
        # MPI admin will distribute jobs to workers
        run_admin()
    else:
        # Workers will cary out work
        run(mpi_rank)

def print_sub_to_check(axisnum):
    axis_params = gen.list_of_axes[axisnum]
    print(find_sub(str(datetime.datetime(*axis_params.obsdate)), axis_params.radar_lat, axis_params.radar_lon, axis_params.radar_el))

if __name__ == '__main__':
    main()
    #run_test()



'''
def save_points(list_of_points, current_num):

    (axis, job_num)

    list_of_job_limits = [523, 1002, 1490]

    jobs_per_axis = [np.ceil(R_M/gen.range_res[0]), ceil(R_M/gen.range_res[1]) + ceil(R_M/gen.range_res[0])]

    jobs_per_axis = 500
    axes = 7
    total_jobs = axes * jobs_per_axis

    current_job = 1030

    current_axis = int(current_job / jobs_per_axis)
    current_job_IN_axis = current_job % jobs_per_axis


    axis = 1

    path = f'{conf.point_path}a{axis}/{current_num:06}.dat'
    with open(path, 'w') as fp:
        for pair in list_of_points:
            fp.write(f'{pair[0][0]},{pair[0][1]} {pair[1][0]},{pair[1][1]}')
im = np.flip(np.array(Image.open(IMAGE_FILE_LOC).convert('L')).T, axis = 1)
mapfunc = intp.interp2d(512, 1024, im)
'''