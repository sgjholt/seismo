#KAPPA AND Q DETERMINATION FOLLOWING ANDERSON AND HUGHES (1989)
#AUTHOR: JAMES HOLT, UNIVERSITY OF LIVERPOOL, j.holt@liverpool.ac.uk

#-------------------------------------------------------------------------------
#IMPORT MODULES 
import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob
from sys import argv
from utils import *
from build_GMdb import calcFAS
#-------------------------------------------------------------------------------
#argv tracker
#argv[0] - script name
#argv[1] - earthquake file name (eg. '20070716101300.kik')
#argv[2] - surface or borehole
#argv[3] - option to choose window for kappa determination
#argv[4] - option to choose geometric mean 

#-------------------------------------------------------------------------------



def main():

    if argv[3] == '-w':
     #whole waveform
    
    elif argv[3] == '-p':
    #p-window

    elif argv[3] == '-s':
    #s-window




def wholeWave():
    if argv[4] == '-G':
        flist = grabRelevent()
        Bin = int(len(flist) / 2) #need to know how many things in list to ...
                               #separate them accordingly (in half).
        flistEW = flist[0 : Bin]
        flistNS = flist[Bin : (2 * Bin)]

        for i in range(0, len(flistEW)):
            FASGeo = processGeo(flistEW[i], flistNS[i])
    else:
        flist = grabRelevent()
        for i in range(0, len(flist)):
            FAS = processNorm(flist[i])



def grabNoise(HostFile):
    dat = np.loadtxt(str(HostFile), skiprows=17)
    head = np.loadtxt(str(HostFile) + '.h')
    dat = preprocessdata(dat, head)
    
    fig, ax = plt.subplots()
    ax.set_title('custom picker for line data')
    line, = ax.plot(dat, 'o', picker=line_picker)
    fig.canvas.mpl_connect('pick_event', onpick2)
    



def grabRelevant():

    path = "/Volumes/J_Holt_HDD/MRes/Modules/Thesis/Data/"

    flist = grab_file_names(str(path+argv[1]), 2) #grabs all of the file names
    
    if argv[2] == 'surface':
        if argv[4] == '-G':
            flistEW = sepBYcomponant(flist, 'EW2')
            flistNS = sepBYcomponant(flist, 'NS2')
            flist = flistEW+flistNS
        elif argv[4] == '-E':
            flist = sepBYcomponant(flist, 'EW2')
        elif argv[4] == '-N':
            flist = sepBYcomponant(flist, 'NS2')
        elif argv[4] == '-U':
            flist = sepBYcomponant(flist, 'UD2')
        else:
            print('Incorrect')
    elif argv[2] == 'borehole':
        if argv[4] == '-G':
            flistEW = sepBYcomponant(flist, 'EW1')
            flistNS = sepBYcomponant(flist, 'NS1')
        elif argv[4] == '-E':
            flist = sepBYcomponant(flist, 'EW1')
        elif argv[4] == '-N':
            flist = sepBYcomponant(flist, 'NS1')
        elif argv[4] == '-U':
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
    FAS = geoMean(FAS1, FAS2)
    return FAS #header just for sampling frequency.

def processNorm(fname):
    #Load the revelent files
    data = np.loadtxt(str(fname), skiprows=17)
    header = np.loadtxt(str(fname)+'.h')
    #assign real units to the seismograms - reshape
    data = preprocessdata(data, header)
    FAS = calcFAS(data, head1)
    return FAS
    


 
    
