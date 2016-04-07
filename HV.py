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


#-------------------------------------------------------------------------------

def main():
    if argv[3] == '-G':
        horzflist = horzWave()
        vertflist = vertWave()
        Bin = int(len(horzflist) / 2) #need to know how many things in list to ...
                                      #separate them accordingly (in half).
        flistEW = horzflist[0 : Bin]
        flistNS = horzflist[Bin : (2 * Bin)]
        n = 0
        for i in range(0, len(horzflist)):
            dathorz = processGeo(flistEW[i], flistNS[i]) 
            datvert = processNorm(vertflist[i]) 
            n += 1
            
            H_V =  dathorz[:,1]/datvert[:,1]
            plt.subplot(Bin,Bin,n)
            plt.loglog(dathorz[:,0],dathorz[:,1])
            plt.loglog(dathorz[:,0],datvert[:,1])
             
            plt.loglog(dathorz[:,0], H_V)
            #plt.subplot(Bin,Bin, i)  
            
            plt.show()











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
    horzDat = geoMean(dat1, dat2)
    FAS = calcFAS(horzDat, head1)
    
    return FAS #header just for sampling frequency.

def processNorm(fname):
    #Load the revelent files
    data = np.loadtxt(str(fname), skiprows=17)
    header = np.loadtxt(str(fname)+'.h')
    #assign real units to the seismograms - reshape
    data = preprocessdata(data, header)
    FAS = calcFAS(data, header)
    return FAS

   
