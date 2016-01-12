#!/usr/bin/python

#SCRIPT TO PROCESS K-NET/KIK-NET DATA
#AUTHOR: JAMES HOLT - UNIVERSITY OF LIVERPOOL
#CREATED: 2015/10/12. LAST EDITED: 2015/10/13

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
script, filename, header = argv
###############################################################################
#COUNT THE NUMBER OF LINES IN INPUT FILE
with open(filename) as f:
        lc = sum(1 for _ in f)
###############################################################################
#SECTION TO LOAD ORIGINAL FILE INTO PYTHON
a = np.loadtxt(filename, skiprows=17)# skips the header information (17 rows)
##############################################################################
#SECTION TO LOAD INTO PYTHON AND DESCRIBE CONTENTS OF INPUT HEADER FILE
# h[0:1:2] = date of quake (YYYY:MM:DD), h[2:4:5] = time of quake (HR:MIN:SEC) ...
# ... h[6] = quake lat (deg), h[7] = quake long, h[8] = depth(km) ...
# ... h[9] = Mw, h[10] = station lat (deg), h[11] = station long (deg) ...
# ... h[12] =  station height (m), h[13:14:15] = record time date (YYYY:MM:DD) ...
# ... h[16:17:18] = record time (HR:MIN:SEC), h[19] = sampling frequency (Hz) ... h[20] = duration time (s)
# ... h[21] = scaling factor numerator (gal), h[22] = scaling ...
# ... factor denominator (gal), h[23] = max acc (gal) 
h = np.loadtxt(header) 

durT = np.arange(0, h[20], (1/h[19]))
##############################################################################
#SECTION TO CREATE DATE-TIME ARRAY FROM INPUT SEISMOGRAM HEADER FILE
# input recording start time from header file (YR, MNTH, DAY, HR, MIN, SEC,...
# ... msec, microsec, nanosec) int to force interger input
t = [int(h[13]), int(h[14]), int(h[15]), int(h[16]), int(h[17]), int(h[18]), 0, 0, 0]
# Must subtract 15 seconds from start time to account for trigger delay ...
# ... stime is a time.mktime() unix timestamp generator
stime = time.mktime(t) - 15
# Endtime is start time + duration in h[20] of header file
etime = stime + h[20]
#create an np array of unix timestamps 
timestamps = np.arange(stime, etime, (1/h[19]))
# convert unix time to python time object
dates=[dt.datetime.fromtimestamp(ts) for ts in timestamps]
# convert python time object to use with matplotlib.dates
datenums=md.date2num(dates)
##############################################################################
#DEFINE SCALING FACTOR NUMERATOR AND DENOMINATOR
#call the numerator and denominator for scf and add 0.0 so it doesnt round to ..
#.. 0 when they are divided have to float them because they're imported as a string which should be useful when using the UTCDateTime function. 
scalingfactornumerator = h[21] + 0.0
scalingfactordenominator = h[22] + 0.0
##############################################################################
#FUNCTION TO PREPROCESS DATA
def preprocessdata(txt, lc, scfN, scfD): #takes the text file data and lc as input
    #calculate reshape factor
    rf = (lc - 17) * 8 #always 8 columns in file
    #transpose the data to reshape in correct order
    trans = txt.transpose()
    #reshape the data into a rf x 1 array
    res = np.reshape(trans, rf, 1)
    #detrend the data (linear least squares by default)
    det = sg.detrend(res)
    #define scaling factor
    scf = scfN / scfD
    #apply scaling factor to have acceleration data
    real = det * scf
    #return as a numpy array divide by 100 to go from gals to g (9.81 m/s^2)
    return real / 100
    

#the function reshapes the data, detrends and applies the scaling factor ...
# in gals, returns array in g's.
acc = preprocessdata(a, lc, scalingfactornumerator, scalingfactordenominator)
##############################################################################
#PLOT THE RESULTING SEISMOGRAMS
plt.subplot(2, 1, 2)
plt.subplots_adjust(bottom=0.2)
plt.xticks( rotation=25 )
ax = plt.gca()
xfmt = md.DateFormatter('%H:%M:%S')
ax.xaxis.set_major_formatter(xfmt)
plt.plot(datenums,acc)
plt.xlim([np.min(datenums), np.max(datenums)])
plt.xlabel('Time [HR:MIN:SEC]')
plt.ylabel('Acc [g]')

plt.subplot(2, 1, 1)
plt.plot(durT, acc)
plt.xlabel('Time [s]')
plt.ylabel('Acc [g]')
plt.show()
##############################################################################
#Calculate the difference in longitude and latitude and give the absolute value # ...(order of subtraction is not really important here). 
dlat = np.absolute((h[10] + 0.0) - (h[6] + 0.0))
dlon = np.absolute((h[11] + 0.0) - (h[7] + 0.0))

#Convert the difference in lat/lon according to the WGS84 reference spheroid ...
# ... to meters
#dlatm = 111132.954 - (559.822 * np.cos(2 * dlat)) + (1.175 * np.cos(4 * dlat))

#maximum acceleration
n = np.max(np.absolute(acc))
print 'PGA %r g' % (n)
print "Delta Lat = %r and Delta Lon = %r " % (dlat, dlon)


