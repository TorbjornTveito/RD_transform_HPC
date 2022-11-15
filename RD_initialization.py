import datetime
from try_spicey import find_sub
from generator_parameters import list_of_axes

'''
entire file is useless, atm...
'''

def find_srps(axisparams, dt = 60):
    obsdate1 = datetime.datetime(*axisparams.obsdate)
    srp1 = find_sub(str(obsdate1), axisparams.radar_lat, axisparams.radar_lon, axisparams.radar_el)
    obsdate2 = obsdate1 + datetime.timedelta(seconds = dt)
    srp2 = find_sub(str(obsdate2), axisparams.radar_lat, axisparams.radar_lon, axisparams.radar_el)
    return(srp1, srp2)

for axis in list_of_axes:
    srp1, srp2 = find_srps(axis)
    print(f"srp1: {srp1}, srp2: {srp2}" )