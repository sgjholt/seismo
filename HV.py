#H/V RATIO SCRIPT FOLLOWING NAKUMURA METHOD (1989)

#-------------------------------------------------------------------------------

#AUTHOR: JAMES HOLT, UNIVERSITY OF LIVERPOOL - j.holt@liverpool.ac.uk
#DATE CREATED: 2016/04/05

#-------------------------------------------------------------------------------
#IMPORT MODULES
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sc
from utils import *
from sys import argv
from glob import glob
import os

#-------------------------------------------------------------------------------

#the main() function will run the script - the user will define what they ...
 #require by passing arguments in using the argv function at the command line.

#-------------------------------------------------------------------------------
#argv tracker
#argv[0] - script name
#argv[1] - earthquake folder name (eg. '20070716101300.kik')
#argv[2] - surface or borehole
#argv[3] - geometric or no '-G' = geometric
#argv[4] - define upper and lower limit of plots (Max 5)
#argv[5] - define upper limit for frequency range of interest
#argv[6] - smoothing 'on' or 'off'
#-------------------------------------------------------------------------------



def main():
    if argv[3] == '-G':
        horzflist = horzWave()
        vertflist = vertWave()
        Bin = int(len(horzflist) / 2) #need to know how many things in list to ...
                                      #separate them accordingly (in half).
        flistEW = horzflist[0 : Bin]
        flistNS = horzflist[Bin : (2 * Bin)]
        names = grabSiteName(flistEW)
        
        if argv[6] == 'on':
            print('Please define a smoothing window \n') 
            smooth_wind = input('> ')
        else:
            smooth_wind = 0

        n = 0
        for i in range(0, int(argv[4])):
            dathorz, freq = processGeo(flistEW[i], flistNS[i]) 
            datvert = processNorm(vertflist[i]) 
            H_V =  dathorz/datvert
            
            n += 1
            
           
            plt.subplot(int(argv[4]),1, n)
            
            plotHV(dathorz, datvert, H_V, argv[6], names[i], smooth_wind, freq)
        plt.show()

        

            
    

def grabSiteName(flist):
    if argv[2] == 'borehole':

        flist = [f.replace('/Volumes/J_Holt_HDD/MRes/Modules/Thesis/Data/' + str(argv[1]), '') for f in flist]        

        flist = [f.replace('/', '') for f in flist] 
    
        flist = [f.replace('.EW1', '') for f in flist] 

    elif argv[2] == 'surface':

        flist = [f.replace('/Volumes/J_Holt_HDD/MRes/Modules/Thesis/Data/' + str(argv[1]), '') for f in flist]        

        flist = [f.replace('/', '') for f in flist] 
    
        flist = [f.replace('.EW2', '') for f in flist] 
    
    return flist


def runningMeanFast(x, N):
    return np.convolve(x, np.ones((N,))/N)[(N-1):]


def plotHV(horzdat, vertdat, HV, flag, name, smooth_wind, freq):
   

    if flag == 'off':
    
        #plt.loglog(freq[freq<=int(argv[5])],horzdat[freq<=int(argv[5])],label='Horz')
        #plt.loglog(freq[freq<=int(argv[5])],vertdat[freq<=int(argv[5])],label='Vert')
        plt.loglog(freq[freq<=int(argv[5])], HV[freq<=int(argv[5])], label='H/V')
        plt.title('HV Ratio for ' + str(name) + ' :Smoothing OFF')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Amplitude [m/s]')
        plt.legend(loc='upper left', shadow=True)
        

    elif flag == 'on':
        

        freq = runningMeanFast(freq, int(smooth_wind))
        horzdat = runningMeanFast(horzdat, int(smooth_wind))
        vertdat = runningMeanFast(vertdat, int(smooth_wind))
        HV = runningMeanFast(HV, int(smooth_wind))
        
        #plt.loglog(freq[freq<=int(argv[5])],horzdat[freq<=int(argv[5])],'-',label='Horz')
        #plt.loglog(freq[freq<=int(argv[5])],vertdat[freq<=int(argv[5])],'-',label='Vert')
        plt.loglog(freq[freq<=int(argv[5])], HV[freq<=int(argv[5])],'-', label='H/V')
        plt.title('Smoothed HV Ratio for ' + str(name) + ' : Smoothing window =' + str(smooth_wind))
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Amplitude [m/s]')
        plt.legend(loc='upper left', shadow=True)
        
        






def horzWave():
    if argv[3] == '-G':
        flist = grabRelevent()
    else:
        flist = grabRelevent()
    return flist


def vertWave():
    path = "/Volumes/J_Holt_HDD/MRes/Modules/Thesis/Data/"

    flist = grab_file_names(str(path+argv[1]), 2) #grabs all of the file names
    
    if argv[2] == 'surface':
        flist = sepBYcomponant(flist, 'UD2') 
    elif argv[2] == 'borehole':
        flist = sepBycomponant(flist, 'UD1')
    return flist
    


def grabRelevent():

    path = "/Volumes/J_Holt_HDD/MRes/Modules/Thesis/Data/"

    flist = grab_file_names(str(path+argv[1]), 2) #grabs all of the file names
    
    if argv[2] == 'surface':
        if argv[3] == '-G':
            flistEW = sepBYcomponant(flist, 'EW2')
            flistNS = sepBYcomponant(flist, 'NS2')
            flist = flistEW+flistNS
        elif argv[3] == '-E':
            flist = sepBYcomponant(flist, 'EW2')
        elif argv[3] == '-N':
            flist = sepBYcomponant(flist, 'NS2')
        elif argv[3] == '-U':
            flist = sepBYcomponant(flist, 'UD2')
        else:
            print('Incorrect')
    elif argv[2] == 'borehole':
        if argv[3] == '-G':
            flistEW = sepBYcomponant(flist, 'EW1')
            flistNS = sepBYcomponant(flist, 'NS1')
        elif argv[3] == '-E':
            flist = sepBYcomponant(flist, 'EW1')
        elif argv[3] == '-N':
            flist = sepBYcomponant(flist, 'NS1')
        elif argv[3] == '-U':
            flist = sepBYcomponant(flist, 'UD1')
        else:
            print('Incorrect')
    return flist


def processGeo(fname1, fname2):
    #Load the relevent files
    dat1 = np.loadtxt(str(fname1), skiprows=17)
    head1 = np.loadtxt(str(fname1) + '.h')
    dat2 = np.loadtxt(str(fname2), skiprows=17)
    head2 = np.loadtxt(str(fname2) + '.h')
    #assign real units to the seismograms - reshape
    dat1 = preprocessdata(dat1, head1)
    dat2 = preprocessdata(dat2, head2)
    #calculate the geometric mean for the two horizontal axes
    FAS1 = calcFAS(dat1, head1)
    FAS2 = calcFAS(dat2, head2)
    FAS = geoMean(FAS1[:,1], FAS2[:,1])
    freq = FAS1[:,0]
    return FAS, freq #header just for sampling frequency.

def processNorm(fname):
    #Load the revelent files
    data = np.loadtxt(str(fname), skiprows=17)
    header = np.loadtxt(str(fname)+'.h')
    #assign real units to the seismograms - reshape
    data = preprocessdata(data, header)
    FAS = calcFAS(data, header)
    return FAS[:,1]


if argv[0] == 'HV.py':
    main() 
