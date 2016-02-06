#!/usr/bin/python

#FUNCTIONS TO PROCESS/PLOT K/KIK-NET DATA
#AUTHOR: JAMES HOLT - UNIVERSITY OF LIVERPOOL
#CREATED: 2015-10-25 19:00:23. LAST EDITED: 2015-10-28 21:16:10  

###############################################################################
# IMPORT REQUIRED MODULES
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from scipy import signal as sg
from scipy import integrate as igt
from scipy import fft
from sys import argv
from obspy import UTCDateTime
import datetime as dt
import time
###############################################################################
# GIVE VARIABLE NAME TO ARGUMENT VALUES
script, filename, header = argv
###############################################################################
#COUNT THE NUMBER OF LINES IN INPUT FILE
def linecounter(fname):
    with open(fname) as f:
        linecount = sum(1 for _ in f)
        name = f.name
    return linecount, name
###############################################################################
#SECTION TO LOAD ORIGINAL FILE INTO PYTHON
def loadknet_filenhead(fname, fhdr):
    a = np.loadtxt(fname, skiprows=17) #load modified kik/k-net data
    b = np.loadtxt(fhdr)# skips the header information (17 rows)
    return a, b
###############################################################################
#FUNCTION TO PREPROCESS DATA - gives back acceleration time series
def preprocessdata(f, hdr): 
    #calculate reshape factor
    rf = len(f) * 8 #always 8 columns in file
    #tranpose, reshape and detrend the data
    det = sg.detrend(np.reshape((f.transpose()), rf, 1))
    #DEFINE SCALING FACTOR & NUMERATOR/DENOMINATOR
    #call the numerator and denominator for scf and add 0.0 so it doesnt ...
    #...round to 0 when they are divided.    
    scfN = float(hdr[21]) 
    scfD = float(hdr[22])  
    scf = scfN / scfD 
    #apply scaling factor to return acceleration values
    real = det * scf
    #return as a numpy array divide by 100 to go from gals to m/s^2
    return real / 100  
##############################################################################
def definestime(tseries, hdr):
    durT = np.arange(0, ((np.size(tseries)) / hdr[19]), (1 / hdr[19])) 
    # this creates and array based on the size of the input time series ...
    #... the samping rate associated with it.        
    t = [int(hdr[13]), int(hdr[14]), int(hdr[15]), int(hdr[16]),
    int(hdr[17]), int(hdr[18]), 0, 0, 0]
    # Must subtract 15 seconds from start time to account for trigger delay ...
    # ... stime is a time.mktime() unix timestamp generator
    stime = time.mktime(t) - 15
    # Endtime is start time + duration in h[20] of header file
    etime = stime + (np.max(durT))
    #create an np array of unix timestamps 
    timestamps = np.arange(stime, etime, (1 / hdr[19]))
    # convert unix time to python time object
    dates = [dt.datetime.fromtimestamp(ts) for ts in timestamps]
    # convert python time object to use with matplotlib.dates
    datenums = md.date2num(dates)
    return durT, datenums
##############################################################################
def plottingacc(tseries, time, dates, station):
    plt.figure()    
    plt.subplot(2, 1, 1)
    plt.plot(time, tseries)
    plt.xlabel('Time [s]')
    plt.ylabel('Acc [m/s/s]')
  
    
    plt.subplot(2, 1, 2)
    plt.subplots_adjust(bottom=0.2)
    plt.xticks( rotation=25 )
    ax = plt.gca()
    xfmt = md.DateFormatter('%H:%M:%S')
    ax.xaxis.set_major_formatter(xfmt)
    plt.plot(dates, tseries)
    plt.xlim([np.min(dates), np.max(dates)])
    plt.xlabel('Time [HRS:MIN:SEC]')
    plt.ylabel('Acc [m/s/s]')
    plt.suptitle('%s' % station)
    
    
    plt.show()
##############################################################################
def int_and_append(tseries, hdr):
    #integrate the time series using scipy's cumulative trapezium method    
    inted_tseries = igt.cumtrapz(tseries, x=None, dx=(1 / hdr[19]))
    #append a 0 to the end to make array size correct again
    inted_tseries = np.append(inted_tseries, 0)
    return inted_tseries
##############################################################################
def plottingaccveldisp(acc, vel, disp, dates, station):
    
    ax1 = plt.subplot(3, 1, 1)
    ax = plt.gca()
    xfmt = md.DateFormatter('%H:%M:%S')
    ax.xaxis.set_major_formatter(xfmt)
    plt.plot(dates, acc, 'k')
    plt.xlim([np.min(dates), np.max(dates)])
    plt.ylabel('Acc [m/s/s]')
    plt.setp(ax1.get_xticklabels(), visible = False)
    

    ax2 = plt.subplot(3, 1, 2, sharex=ax1)
    plt.xticks( rotation=25 )
    plt.plot(dates, vel, 'r')
    plt.xlim([np.min(dates), np.max(dates)])
    plt.ylabel('Vel [m/s]')
    plt.setp(ax2.get_xticklabels(), visible = False)

    ax3 = plt.subplot(3, 1, 3, sharex=ax1)
    plt.plot(dates, disp, 'b')
    plt.xlabel('Time [HRS:MIN:SEC]')
    plt.xlim([np.min(dates), np.max(dates)])
    plt.ylabel('Disp [m]')
    plt.suptitle('%s' % station)
    plt.setp(ax3.get_xticklabels())
    plt.subplots_adjust(bottom=0.2)
    plt.xticks( rotation=25 )
    
    plt.show()
##############################################################################
def frequencyplots(tseries, hdr, station):
    #calculate freq range for data
    Fs = hdr[19] 
    n = len(tseries) #length of the signal
    k = np.arange(n) #build an array based on length of signal
    T = n/Fs #sampling time delta
    frq = k/T #two sides frequency range
    freq = frq[range(n/2)] #one side frequency range
    #fft computing and normalisation with hamming window applied to ...
    #spectral leakage 
    Y = np.fft.fft(tseries) / n 
    Z = Y[range(n/2)] #one side amplitude range
    
  
    plt.subplot(3,1,1)
    plt.plot(freq, abs(Z), 'g-')
    plt.xlabel('Freq [Hz]')
    plt.ylabel('|Spectral Acc [m/s]|')
    plt.grid()
    
    plt.subplot(3,1,2)
    plt.semilogy(freq, abs(Z), 'r-')
    plt.xlabel('Freq [Hz]')
    plt.ylabel('Log10(|Spectral Acc [m/s]|)')
    plt.grid()

    plt.subplot(3,1,3)
    plt.loglog(freq, abs(Z), 'k-')
    plt.xlabel('Log10(Freq [Hz])')
    plt.ylabel('Log10(|Spectral Acc [m/s]|)')
    plt.grid()  
    plt.axis('tight')
    plt.suptitle('Amp Spectra/Fourier Amp Spectra for %s' % station)   
    plt.tight_layout()
    plt.show() 

############################################################################

f, hdr = loadknet_filenhead(filename, header)

lc, name = linecounter(filename)

acc = preprocessdata(f, lc, hdr)

t, d = definestime(acc, hdr)

plottingacc(acc, t, d, name)

vel = int_and_append(acc, hdr)

disp = int_and_append(vel, hdr)

plottingaccveldisp(acc, vel, disp, d, name)

frequencyplots(acc, hdr, name)


