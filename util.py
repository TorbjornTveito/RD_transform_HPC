import pickle
import os
import numpy as np
import config as conf


def cartesian2lunar(c):
    north=np.array([0.0,0.0,1.0])
    x=np.array([1.0,0.0,0.0])
    y=np.array([0.0,1.0,0.0])
    lat=90.0-180.0*np.arccos(np.dot(c,north)/(np.sqrt(np.dot(c,c))))/np.pi

    eq=np.dot(c,x)*x + np.dot(c,y)*y
    lon=180.0*np.arccos(np.dot(x,eq)/(np.sqrt(np.dot(eq,eq)*np.sqrt(np.dot(x,x)))))/np.pi
    if np.dot(c,y) < 0.0:
        lon=360.0-lon
    return(lat,lon)

def lunar2cartesian(lat,lon):
    return(np.array([np.cos(np.pi*lat/180.0)*np.cos(np.pi*lon/180.0), (np.cos(np.pi*lat/180.0)*np.sin(np.pi*lon/180.0)), (np.sin(np.pi*lat/180.0))]))
