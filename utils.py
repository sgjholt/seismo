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

    if flag == 2:
        f1 = glob(str(path_to_folder) + '*.EW1')

        f2 = glob(str(path_to_folder) + '*.NS1')

        f3 = glob(str(path_to_folder) + '*.UD1')

        f4 = glob(str(path_to_folder) + '*.EW2')
    
        f5 = glob(str(path_to_folder) + '*.NS2')

        f6 = glob(str(path_to_folder) + '*.UD2')

        joinedFileList =  f1 + f2 + f3 + f4 + f5 + f6
    
    return joinedFileList

def preprocessdata(f, hdr): 
    #calculate reshape factor
    rf = len(f) * 8 #always 8 columns in file
    #transpose, reshape and detrend the data
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
    

    
    geoMean = (file1*file2)**(1/2)

    return geoMean

def sepBYcomponant(fnames, flag):

    Bin = int(len(fnames)/6) 

    if flag == 'EW1':
        fnames = fnames[0:Bin]
    elif flag == 'NS1':
        fnames = fnames[(Bin+1):(2*Bin)]
    elif flag == 'UD1':
        fnames = fnames[(2*Bin+1):(3*Bin)]
    elif flag == 'EW2':
        fnames = fnames[(3*Bin+1):(4*Bin)]
    elif flag == 'NS2':
        fnames = fnames[(4*Bin+1):(5*Bin)]
    elif flag == 'UD2':
        fnames = fnames[(5*Bin+1):(6*Bin)]
    else:
        print("Incorrect Option.")
    
    return fnames

def calcFAS(tseries, h):
    Fs = h[19]
    n = len(tseries) #length of the signal
    k = np.arange(n) #build an array based on length of signal
    T = n/Fs #sampling time delta
    frq = k/T #two sides frequency range
    freq = frq[range(int(n/2))] #one side frequency range
    #fft computing and normalisation with hamming window applied to ...
    #spectral leakage ls

    Y = np.fft.fft(tseries) / n 
    Z = Y[range(int(n/2))] #one side amplitude range    
    FAS = [freq, abs(Z)]
    return np.transpose(FAS)

#def line_picker(line, mouseevent):
 #   """
  #  find the points within a certain distance from the mouseclick in
   # data coords and attach some extra attributes, pickx and picky
    #which are the data points that were picked
   # """
    #if mouseevent.xdata is None:
    #    return False, dict()
    #xdata = line.get_xdata()
    #ydata = line.get_ydata()
    #maxd = 0.05
    #d = np.sqrt((xdata - mouseevent.xdata)**2. + (ydata - mouseevent.ydata)**2.)

    #ind = np.nonzero(np.less_equal(d, maxd))
    #if len(ind):
       # pickx = np.take(xdata, ind)
      #  picky = np.take(ydata, ind)
     #   props = dict(ind=ind, pickx=pickx, picky=picky)
        #return True, props
    #else:
       # return False, dict()

#def onpick2(event):
 #   print('onpick2 line:', event.pickx, event.picky)












    
