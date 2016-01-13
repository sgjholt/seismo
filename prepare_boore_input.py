
import numpy as np
from glob import glob
from sys import argv
from process_kikhdr import W_Cempfaultparams, hyper_epiD
import os.path


#path_to_folder = argv[1]

#filenames1 = glob('/Users/jamesholt/Documents/MRes/20001006133000/20001006133000.kik/*.EW1.h')

#filenames2 = glob('/Users/jamesholt/Documents/MRes/20001006133000/20001006133000.kik/*.NS1.h')

#filenames3 = glob('/Users/jamesholt/Documents/MRes/20001006133000/20001006133000.kik/*.UD1.h')
#---------------

# argv in in order of argv[0] = python script name, argv[1] = full folder path, 
#... argv[2] is dip in degrees, argv[3] is the magnitude (Mw) and argv[4] 
#... is the output file name.

def main():
    path_to_folder = str(argv[1])   
    #e.g. '/Users/jamesholt/Documents/MRes/20001006133000/20001006133000.kik/'
    
    
    fileList = grab_file_names(path_to_folder)
    for fname in fileList:
        make_db_files(fname)
    
    reshapeFile()
    
    saveEqParams(path_to_folder)
        
    write_fname_to_file(fileList)
    
    
    
        
      


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
    with open (str(argv[1]) + str(argv[4]), 'ab') as f:          
            
        h = np.loadtxt(fname)
        EpiD, HypD = hyper_epiD(h)
        
        dbfilearray = np.array([h[10], h[11], h[12], h[23]/100, EpiD/1000, HypD/1000])
        np.savetxt(f, dbfilearray)
    
    

def reshapeFile():
    
    dat = np.loadtxt(str(argv[1]) + str(argv[4]))
    np.savetxt(str(argv[1]) + str(argv[4]), dat.reshape(len(dat) / 6, 6))

def saveEqParams(path_to_folder):
    
    with open(path_to_folder + 'quake_params.txt', 'wb') as f:
        
    
        dip, mag = [float(argv[2]), float(argv[3])]
        
        SL, DDW, DDWA, RSA, SW, SWA  = W_Cempfaultparams(mag, dip)
    
        List = grab_file_names(path_to_folder)
        h = np.loadtxt(List[0])
        eqparams = np.array([h[6], h[7], h[8], dip, mag, SL, DDW, DDWA, RSA, SW, SWA])  
            
        header = "Quake lat(deg), Quake lon(deg), Quake Depth(km), Quake Dip (deg), Quake Mag (Mw) Surface Length(km), Downdip Width(km) (with formula), Downdip W(km) (with area), Area(km^2), S width w/out area(km), S width with area(km)"         
        np.savetxt(f, eqparams, fmt='%10.5f', header=header)


main()

        
    

#def check_if_exist(path_to_file, filename):
    #if os.path.isfile(path_to_file + filename) == False:
        #break
    #else: 
        #print('Warning file already exists, please remove or rename.')
        

          
   

            
            











