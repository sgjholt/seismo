#!/usr/bin/python

#SCRIPT TO PROCESS K-NET/KIK-NET HEADER FILE DATA
#AUTHOR: JAMES HOLT - UNIVERSITY OF LIVERPOOL
#CREATED: 2015/10/19. LAST EDITED: 2015-11-02 21:10:26 

###############################################################################
# IMPORT REQUIRED MODULES
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from scipy import signal as sg
from sys import argv
from obspy import UTCDateTime
import datetime as dt
import time
###############################################################################
# GIVE VARIABLE NAME TO ARGUMENT VALUES
script, header, fdip, Mag = argv
###############################################################################
#COUNT THE NUMBER OF LINES IN INPUT FILE
with open(header) as f:
        fname = f.name
###############################################################################
#SECTION TO LOAD INTO PYTHON AND DESCRIBE CONTENTS OF INPUT HEADER FILE
# h[0:1:2] = date of quake (YYYY:MM:DD), h[2:4:5] = time of quake (HR:MIN:SEC) ...
# ... h[6] = quake lat (deg), h[7] = quake long, h[8] = depth(km) ...
# ... h[9] = Mw, h[10] = station lat (deg), h[11] = station long (deg) ...
# ... h[12] =  station height (m), h[13:14:15] = record time date (YYYY:MM:DD) ...
# ... h[16:17:18] = record time (HR:MIN:SEC), h[19] = sampling frequency (Hz) ... 
# h[20] = duration time (s) h[21] = scaling factor numerator (gal), 
# h[22] = scaling factor denominator (gal), h[23] = max acc (gal) 

h = np.loadtxt(header) 

##############################################################################
#Calculate the difference in longitude and latitude and give the absolute value 
# ...(order of subtraction is not really important here because distance is ...
# so short.). 
def hyper_epiD(hdr):    
    
    dlat = np.absolute((hdr[10] + 0.0) - (hdr[6] + 0.0))
    dlon = np.absolute((hdr[11] + 0.0) - (hdr[7] + 0.0))
    #Using WGS84 Spheroid Determin what 1 degree is in lat/lon in meters at ...
    # ... at seismic network location (should not change much)
    stnlat = hdr[10] * (np.pi / 180)
    onedlat = 111132.92 - (559.82 * np.cos(2 * stnlat)) + (1.175 * np.cos(4 * stnlat)) - (0.0023 * np.cos(6 * stnlat))
    onedlon = (111412.84 * np.cos(stnlat)) - (93.5 * np.cos(3 * stnlat)) - (0.118 * np.cos(5 * stnlat))
    # multipy dlat/dlon by 1 degree lat/lon to get distance in m
    latm = dlat * onedlat
    lonm = dlon * onedlon
    # calculate point source 2d epicentral distance in m
    EpiD = np.sqrt((latm * latm) + (lonm * lonm)) 
    #calculate point source 2d hypocentral distance in m
    dZ =  (hdr[8] * 1000) + hdr[12]
    HypD = np.sqrt((latm * latm) + (lonm * lonm) + (dZ * dZ)) 
    
    return EpiD, HypD


#Wells and Coppersmith (1994) emprical relationship to determine fault plane ...
#... parameters.


def W_Cempfaultparams(Mag, dip):
    #surface length (KM)
    SRL = 10**((-3.22+0.69*Mag)) 
    #downdip width (KM)
    RW = 10**((-1.01+0.32*Mag))
    #rupture surface area (KM^2)
    RA = 10**((-3.49+0.91*Mag))
    #Downdip width calculated with area / surface length (KM)
    RWA = RA / SRL
    #surface width (KM) Calculated using downdip width derived from... 
    #...rupture area / surface length
    SWA = np.cos((dip * (np.pi / 180))) * RWA
    #surface width (KM)  calculated with direct rupture width
    SW = np.cos((dip * (np.pi / 180))) * RW

    
    return SRL, RW, RWA, RA, SW, SWA


    

#maximum acceleration (from header file / 100 to get to units of m/s/s)
amax = h[23] / 100
#epicentral and hypocentral distances
epid, hypd = hyper_epiD(h)
#empirical fault parameters
srl, rw, rwa, ra, sw, swa = W_Cempfaultparams(float(Mag), float(fdip))


print "%r %r %r %r %r %r" % (h[10], h[11], h[12], amax, (epid / 1000), (hypd / 1000))

b = np.array([h[6], h[7], h[8], srl, rw, rwa, ra, sw, swa])
np.savetxt('quakeparams.txt', b, fmt='%10.5f', header="Quake lat(deg), Quake lon(deg), Quake Depth(km), Surface Length(km), Downdip Width(km) (with formula), Downdip W(km) (with area), Area(km^2), S width w/out area(km), S width with area(km)")




