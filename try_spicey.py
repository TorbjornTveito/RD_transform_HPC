import spiceypy as spc
import math
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time, TimeDelta
from datetime import datetime, timedelta

spc.furnsh('./data/generator_input/kernels.tm')

def vector_length(vec, dim = 3):
    partsum = 0
    for i in range(dim):
        partsum += vec[i]**2
    final = np.sqrt(partsum)
    return(final)

def rev(angle):
    if angle < 0:
        return (angle % 360) - 360
    else:
        return angle % 360

def vector_rotator(vec, frame1, frame2, time):

    vec_pointer = vec / vector_length(vec, len(vec))
    rotation_matrix = spc.pxform(frame1, frame2, time)
    rotated_vector = spc.mxv(rotation_matrix, vec_pointer)
    return(-rotated_vector)

def find_sub(datestring, lat = 69.340, lon = 20.313, el = 0.1): #returns SRP in form of (lon, lat) at given time
    Obstime = spc.str2et(datestring)
    r_moon = 1737.4
    lon = lon *spc.rpd()
    lat = lat *spc.rpd()

    earth_shape_stuff = spc.bodvrd('earth', 'RADII', 3)

    f = (earth_shape_stuff[1][0] - earth_shape_stuff[1][2])/earth_shape_stuff[1][0]

    coords = spc.georec(lon, lat, el, earth_shape_stuff[1][0], f)

    obs_to_moon_center, LTmoon = spc.spkcpo('moon', Obstime, 'ITRF93', 'observer', 'NONE', coords, 'EARTH', 'ITRF93')

    obs_to_MC_xyz = [obs_to_moon_center[0], obs_to_moon_center[1], obs_to_moon_center[2]]

    moon_to_obs_pos = vector_rotator(obs_to_MC_xyz, 'ITRF93', 'MOON_ME', Obstime)

    moon_to_obs_pos *= r_moon

    rad, lon, lat = spc.reclat(moon_to_obs_pos)

    lon = lon*180/np.pi
    if lon > 360:
        lon -= 360
    lat = lat*180/np.pi
    return(lon, lat)


def sub_diff(datestring, timedelta = 60): #returns vector from one SRP to another with time difference timedelta (in seconds)
    Obstime = spc.str2et(datestring)
    r_moon = 1737.4


    e3d_lat = 69.340 * spc.rpd()
    e3d_lon = 20.313 * spc.rpd()

    earth_shape_stuff = spc.bodvrd('earth', 'RADII', 3)

    f = (earth_shape_stuff[1][0] - earth_shape_stuff[1][2])/earth_shape_stuff[1][0]

    e3d_coords = spc.georec(e3d_lon, e3d_lat, 0.1, earth_shape_stuff[1][0], f)

    obs_to_moon_center, LTmoon = spc.spkcpo('moon', Obstime, 'ITRF93', 'observer', 'NONE', e3d_coords, 'EARTH', 'ITRF93')

    obs_to_MC_xyz = [obs_to_moon_center[0], obs_to_moon_center[1], obs_to_moon_center[2]]

    obs2_to_moon_center, LTmoon2 = spc.spkcpo('moon', Obstime + timedelta, 'ITRF93', 'observer', 'NONE', e3d_coords, 'EARTH', 'ITRF93')

    obs2_to_MC_xyz = [obs2_to_moon_center[0], obs2_to_moon_center[1], obs2_to_moon_center[2]]

    moon_to_obs_pos = vector_rotator(obs_to_MC_xyz, 'ITRF93', 'MOON_ME', Obstime)

    moon_to_obs_pos *= r_moon

    moon_to_obs_pos2 = vector_rotator(obs2_to_MC_xyz, 'ITRF93', 'MOON_ME', Obstime + timedelta)

    moon_to_obs_pos2 *= r_moon


    rad, lon, lat = spc.reclat(moon_to_obs_pos)
    rad2, lon2, lat2 = spc.reclat(moon_to_obs_pos2)

    diffvec = (moon_to_obs_pos2 - moon_to_obs_pos)
    diffvec = np.array(diffvec)

    xprod = np.cross(diffvec, moon_to_obs_pos)

    lon = lon*180/np.pi +360
    if lon > 360:
        lon -= 360
    lat = lat*180/np.pi
    lon2 = lon2*180/np.pi +360
    if lon2 > 360:
        lon2 -= 360
    lat2 = lat2*180/np.pi



def find_np(Obstime):
    e3d_lat = 69.340 * spc.rpd()
    e3d_lon = 20.313 * spc.rpd()

    Moon_radii = spc.bodvrd('moon', 'RADII', 3)

    f_moon = (Moon_radii[1][0] - Moon_radii[1][2])/Moon_radii[1][0]


    NP_coords = spc.georec(90 * spc.rpd(), 0, 0.1, Moon_radii[1][0], f_moon)



    NP_to_earth_center, LTmoon = spc.spkcpo('EARTH', Obstime, 'ITRF93', 'observer', 'NONE', NP_coords, 'moon', 'MOON_ME')

    Earth_radii = spc.bodvrd('earth', 'RADII', 3)

    f_earth = (Earth_radii[1][0] - Earth_radii[1][2])/Earth_radii[1][0]

    e3d_coords = spc.georec(e3d_lon, e3d_lat, 0.1, Earth_radii[1][0], f_earth)

    NP_to_e3d = [NP_to_earth_center[0] + e3d_coords[0], NP_to_earth_center[1] + e3d_coords[1], NP_to_earth_center[2] + e3d_coords[2]]

def find_dist(date, lat, lon, el):
    r_moon = 1737.4

    lat = 69.58 * spc.rpd()
    lon = 19.23 * spc.rpd()

    earth_shape_stuff = spc.bodvrd('earth', 'RADII', 3)

    f_earth = (earth_shape_stuff[1][0] - earth_shape_stuff[1][2])/earth_shape_stuff[1][0] #flattening coefficient for earth

    coords = spc.georec(lon, lat, el, earth_shape_stuff[1][0], f_earth)

    Moon_radii = spc.bodvrd('moon', 'RADII', 3)
    f_moon = (Moon_radii[1][0] - Moon_radii[1][2])/Moon_radii[1][0]

    sub = find_sub(date,lat, lon, el)
    obstime = spc.str2et(date)

    #print("sub deg: ",sub)
    sub = [sub[0]*spc.rpd(), sub[1]*spc.rpd()]
    #print("sub rad: ", sub)
    sub_coords = spc.georec(sub[0], sub[1], 0, Moon_radii[1][0], f_moon)
    #print("sub_coords: ",sub_coords)
    sub_to_earth, ltmoon = spc.spkcpo('EARTH', obstime, 'ITRF93', 'observer', 'NONE', sub_coords, 'moon', 'MOON_ME')
    sub_to_e3d = [sub_to_earth[0] + coords[0], sub_to_earth[1] + coords[1], sub_to_earth[2] + coords[2]]
    #print("sub_to_e3d: ", sub_to_e3d)
    obs_to_moon_center, LTmoon = spc.spkcpo('moon', obstime, 'ITRF93', 'observer', 'NONE', coords, 'EARTH', 'ITRF93')
    #print("obs_to_moon_center: ", vector_length(obs_to_moon_center))
    #print("diff: ", vector_length(obs_to_moon_center)-vector_length(sub_to_e3d))
    #print("to sub: ", vector_length(sub_to_e3d))
    return(sub_to_e3d)

#394606.0248616907
if __name__ == "__main__":
    date1 = datetime(2022, 2, 14, 18, 18, 16, 819485)
    ET = spc.datetime2et(date1)

    r_moon = 1737.4
    print(str(date1))

    e3d_lat = 69.58 * spc.rpd()
    e3d_lon = 19.23 * spc.rpd()

    lat = 69.58
    lon = 19.23

    earth_shape_stuff = spc.bodvrd('earth', 'RADII', 3)

    r_e = earth_shape_stuff[1][0]

    f = (earth_shape_stuff[1][0] - earth_shape_stuff[1][2])/earth_shape_stuff[1][0] #flattening coefficient for earth

    find_dist(date1)