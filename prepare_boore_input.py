
import numpy as np
from glob import glob
from sys import argv
from process_kikhdr import W_Cempfaultparams, hyper_epiD
import os
import stat

# argv in in order of argv[0] = python script name, argv[1] = full folder path, 
#... argv[2] is dip in degrees, argv[3] is the magnitude (Mw), argv[4] is the 
#...strike of the fault and argv[5] is the output file name.
#argv[1] e.g. = '/Users/jamesholt/Documents/MRes/20001006133000/20001006133000.kik/'

def main():
    fileList = grab_file_names(str(argv[1]))
    
    for fname in fileList:
        make_db_files(fname)
    
    reshapeFile()
    
    saveEqParams(path_to_folder)
        
    write_fname_to_file(fileList)
    
    booreFileMaker()
    
    
    
        
def grab_file_names(path_to_folder):
    """ This function uses glob to search for file names using the path string 
        to the folder. Glob then searches for the file extension below and 
        stores the filenames in a list. The function retuns a list for each
        of the glob searches separately. """ 
    f1 = glob(path_to_folder + '*.EW1.h')

    f2 = glob(path_to_folder + '*.NS1.h')

    f3 = glob(path_to_folder + '*.UD1.h')

    f4 = glob(path_to_folder + '*.EW2.h')
    
    f5 = glob(path_to_folder + '*.NS2.h')

    f6 = glob(path_to_folder + '*.UD2.h')

    joinedFileList =  f1 + f2 + f3 + f4 + f5 + f6

    return joinedFileList

    

def write_fname_to_file(fname):
    """ This function writes all of the file paths and names to a single file. """
    with open('filenames.txt', 'w') as fnames:   
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
    
    with open(path_to_folder + 'quake_params.txt', 'wb') as f:
        
    
        dip, mag, strike = [float(argv[2]), float(argv[3]), float(argv[4])]
        
        SL, DDW, DDWA, RSA, SW, SWA  = W_Cempfaultparams(mag, dip)
    
        List = grab_file_names(path_to_folder)
        h = np.loadtxt(List[0])
        eqparams = np.array([h[6], h[7], h[8], dip, mag, strike, SL, DDW, DDWA, RSA, SW, SWA])  
            
        header = "Quake lat(deg), Quake lon(deg), Quake Depth(km), Quake Dip (deg), Quake Mag (Mw), Quake Strike (deg), Surface Length(km), Downdip Width(km) (with formula), Downdip W(km) (with area), Area(km^2), S width w/out area(km), S width with area(km)"         
        np.savetxt(f, eqparams, fmt='%10.5f', header=header)


def booreFileMaker():
    fname = str(argv[1]) + 'DIST_3D.CTL'
    
    writeString1 = "!Control file for program Dist_3D\n!Name of output file:\ndist_3d.out\n!Minimum Depth for Campbell:\n3.0\n!Number of Fault Planes:\n1\n!Fault orientation (ref lat,long,elev(m),depth, strike, dip, s1, s2, w1, w2):\n"


    
    writeString2 = "!Lat,long, elev (m) of station, character string for output (can be blank):\n" 

    with open (fname, 'w') as f:
        f.write(writeString1)
  
    with open (fname, 'ab') as f:
        F = np.loadtxt(str(argv[1]) + 'quake_params.txt', skiprows=1)
        boorelist = np.array([F[0], F[1], 0, F[2], F[3], F[5], F[6]/2, F[6]/2, F[11]/2, F[11]/2])
        boorelist = boorelist.reshape(1, 10)
        np.savetxt(f, boorelist)

    with open (fname, 'a') as f:
        f.write(writeString2)  
    
    with open (fname, 'ab') as f:
        File = np.loadtxt(str(argv[1]) + str(argv[5]))
        np.savetxt(f, File[:, [1,2,3]])
   
    with open (fname, 'a') as f:
        f.write('stop')
    
    make_executable(fname)   




def make_executable(path):
    mode = os.stat(path).st_mode
    mode |= (mode & 292) >> 2
    os.chmod(path, mode)



main()

        
    

#def check_if_exist(path_to_file, filename):
    #if os.path.isfile(path_to_file + filename) == False:
        #break
    #else: 
        #print('Warning file already exists, please remove or rename.')
        

          
   

            
            











