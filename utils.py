import numpy as np
import os
from glob import glob
import scipy.signal as sg


def chkEmptyfile(fname):
    for filename in fname:
        f = np.loadtxt(filename)

        if f.size != 0:
            continue
        

        elif f.size == 0:
            print('File {} is empty, please check the source file.'.format(
            filename))
    print('Check is complete.')   

def grab_file_names(path_to_folder, flag):
    """ This function uses glob to search for file names using the path string 
        to the folder. Glob then searches for the file extension below and 
        stores the filenames in a list. The function retuns a list for each
        of the glob searches separately. This function returns all lists added 
        together. """ 
    
    
    if flag == 0:
        f1 = glob(str(path_to_folder) + '*.EW1.h')

        f2 = glob(str(path_to_folder) + '*.NS1.h')

        f3 = glob(str(path_to_folder) + '*.UD1.h')

        f4 = glob(str(path_to_folder) + '*.EW2.h')
    
        f5 = glob(str(path_to_folder) + '*.NS2.h')

        f6 = glob(str(path_to_folder) + '*.UD2.h')

        joinedFileList =  f1 + f2 + f3 + f4 + f5 + f6

    if flag == 1:
        f = glob(str(path_to_folder) + '*.*[0-9]')
        
        joinedFileList = f
    
    return joinedFileList

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

  
def geoMean(file1, file2):
    """Calculates the geometric mean of two seismogram componants (e.g. N-S and
    E-W componants. The files must already be loaded into python using np.load.
    USAGE: geomean = geoMean(EW, NS). """
    file1min = abs(min(file1))
    file2min = abs(min(file2))

    if file1min > file2min:
        Min = file1min
    else:
        Min = file2min
    
    file1 = file1 + (Min + 1)
    file2 = file2 + (Min + 1)

    
    geoMean = (file1*file2)**(1/2)
    geoMean = sg.detrend(geoMean)

    return geoMean
    
