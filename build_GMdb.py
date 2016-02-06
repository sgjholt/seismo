
import numpy as np
from glob import glob
from sys import argv
from process_kikhdr import W_Cempfaultparams, hyper_epiD
import os
import stat
import subprocess
from plot_waveforms import preprocessdata

# argv in in order of argv[0] = python script name, argv[1] = full folder path, 
#... argv[2] is dip in degrees, argv[3] is the magnitude (Mw), argv[4] is the 
#...strike of the fault and argv[5] is the output file name.
#argv[1] e.g. = '/Users/jamesholt/Documents/MRes/20001006133000/20001006133000.kik/'

def main():
    #call all of the filenames that are needed to build databases
    headerfiles = grab_file_names(str(argv[1]), 0)#to build GM database
    files = grab_file_names(str(argv[1]), 1)#to produce FAS for each seismogram
    
    subprocess.call(['rm', str(argv[1]) + str(argv[5])])    
    for fname in headerfiles:
        make_db_files(fname)
    reshapeFile()
    
    saveEqParams(str(argv[1]))
        
    write_fname_to_file(fileList)
    
    booreFileMaker()
    
    runBoore()
    
    mergedatabase()
    
    x = 0
    for i in range(0, len(files)):
        x += 1
        print('Calculating FAS {} out of {}'.format(x, len(files)))
        FASmaker(files[i], headerfiles[i])
    
    
        
def grab_file_names(path_to_folder, flag):
    """ This function uses glob to search for file names using the path string 
        to the folder. Glob then searches for the file extension below and 
        stores the filenames in a list. The function retuns a list for each
        of the glob searches separately. This function returns all lists added 
        together. """ 
    
    
    if flag == 0:
    
    
        f1 = glob(path_to_folder + '*.EW1.h')

        f2 = glob(path_to_folder + '*.NS1.h')

        f3 = glob(path_to_folder + '*.UD1.h')

        f4 = glob(path_to_folder + '*.EW2.h')
    
        f5 = glob(path_to_folder + '*.NS2.h')

        f6 = glob(path_to_folder + '*.UD2.h')

        joinedFileList =  f1 + f2 + f3 + f4 + f5 + f6

    

    if flag == 1:
        f = glob(path_to_folder + '*.*[0-9]')
        joinedFileList = f
    
    return joinedFileList

def write_fname_to_file(fname):
    """ This function writes all of the file paths and names to a single file. """
    subprocess.call(['rm', str(argv[1]) + 'filenames.txt'])
    with open(str(argv[1]) + 'filenames.txt', 'w') as fnames:   
        for f in fname:
            fnames.write(f + '\n')
       

def make_db_files(fname):
    """ This function will create a database from the kik-net headers by 
        extracting stnlat h[10], stnlon h[11], stnheight h[12], max acceleration """
      
    with open (str(argv[1]) + str(argv[5]), 'ab') as f:          
            
        h = np.loadtxt(fname)
        EpiD, HypD = hyper_epiD(h)
        
        dbfilearray = np.array([h[10], h[11], h[12], h[23]/100, EpiD/1000, HypD/1000])
        np.savetxt(f, dbfilearray)
    
    

def reshapeFile():
    """ This function reshapes the prelim database data. Because it won't let me 
        reshape the data in the function that makes it without messing that
        file up. Why? Reasons, thats why. """
    dat = np.loadtxt(str(argv[1]) + str(argv[5]))
    np.savetxt(str(argv[1]) + str(argv[5]), dat.reshape(len(dat) / 6, 6))

def saveEqParams(path_to_folder):
    subprocess.call(['rm', path_to_folder + 'quake_params.txt'])
    with open(path_to_folder + 'quake_params.txt', 'wb') as f:
        
    
        dip, mag, strike = [float(argv[2]), float(argv[3]), float(argv[4])]
        
        SL, DDW, DDWA, RSA, SW, SWA  = W_Cempfaultparams(mag, dip)
    
        List = grab_file_names(path_to_folder)
        h = np.loadtxt(List[0])
        eqparams = np.array([h[6], h[7], h[8], dip, mag, strike, SL, DDW, DDWA, RSA, SW, SWA])  
            
        header = "Quake lat(deg), Quake lon(deg), Quake Depth(km), Quake Dip (deg), Quake Mag (Mw), Quake Strike (deg), Surface Length(km), Downdip Width(km) (with formula), Downdip W(km) (with area), Area(km^2), S width w/out area(km), S width with area(km)"         
        np.savetxt(f, eqparams, fmt='%10.5f', header=header)


def booreFileMaker():
    """booreFileMaker creates a control file to input in boores DIST_3D code. """

    #Define the strings to be used for removing, opening and writing to files.
    fname = '/Users/jamesholt/seismograms/DIST_3D.CTL'
    writeString1 = "!Control file for program Dist_3D\n!Name of output file:\ndist_3d.out\n!Minimum Depth for Campbell:\n3.0\n!Number of Fault Planes:\n1\n!Fault orientation (ref lat,long,elev(m),depth, strike, dip, s1, s2, w1, w2):\n"
    writeString2 = "!Lat,long, elev (m) of station, character string for output (can be blank):\n" 
    #Remove the existing DIST_3D.CTL file.
    subprocess.call(['rm', str(fname)])
    #Begin making the DIST_3D.CTL file
    with open (fname, 'w') as f:
        f.write(writeString1)
    #Add the line which will define fault plane
    #Meaning in order of list =  ref lat, lon, elev, depth, strike, dip, ...
    #length of surface rupture left of reference lat, same but right, ...
    #width of fault projected to the surface to left of reference, same but ...
    #right.
    
    with open (fname, 'ab') as f:
        F = np.loadtxt(str(argv[1]) + 'quake_params.txt', skiprows=1)
        boorelist = np.array(
        [F[0], F[1], 0, F[2], F[5], F[3], F[6]/2, F[6]/2, F[11]/2, F[11]/2])
        boorelist = boorelist.reshape(1, 10)
        np.savetxt(f, boorelist, fmt='%10.5f')

    with open (fname, 'a') as f:
        f.write(writeString2)  
    
    with open (fname, 'ab') as f:
        File = np.loadtxt(str(argv[1]) + str(argv[5]))
        np.savetxt(f, File[:, [0,1,2]], fmt='%10.5f')
   
    with open (fname, 'a') as f:
        f.write('stop')
    
    make_executable(fname)   


def make_executable(path):
    mode = os.stat(path).st_mode
    mode |= (mode & 292) >> 2
    os.chmod(path, mode)

def runBoore():
    path_to_exe = "/Users/jamesholt/Documents/dist_programs/DIST_3D.EXE"    
    subprocess.call(["wine", path_to_exe])

def loadOutFile():
    path_to_outfile = "/Users/jamesholt/seismograms/dist_3d.out"
    outfileDat = np.loadtxt(path_to_outfile, skiprows=13)
    booredists = outfileDat[:, [4, 5, 6, 7]]
    return booredists

def mergedatabase():
    subprocess.call = (['rm', str(argv[1]) + str(argv[5]) + '_complete.txt'])
    with open (str(argv[1]) + str(argv[5]) + '_complete.txt', 'wb') as f:
        maindb = np.loadtxt(str(argv[1]) + str(argv[5]))
        boordb = loadOutFile()
        whole = np.concatenate((maindb, boordb), axis = 1)
        header = "Stn Lat (deg), Stn Lon (deg), Stn Height (m), max acc (m/s/s), epicentral dist (km), hypocentral dist (km), Rrup (km), Rcmpbl (km), Rjb (km), AZjb (deg)"
        
        np.savetxt(f, whole, fmt ='%10.5f', header = header)
    
def FASmaker(fname, hdr):
    
    f = np.loadtxt(fname, skiprows=17)
    h = np.loadtxt(hdr)

    tseries = preprocessdata(f, h)
    FAS = calcFAS(tseries, h)
    F =  (str(fname) + '.fas')
    np.savetxt(F, FAS)



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




main()

        
    

#def check_if_exist(path_to_file, filename):
    #if os.path.isfile(path_to_file + filename) == False:
        #break
    #else: 
        #print('Warning file already exists, please remove or rename.')
        

          
   

            
            











