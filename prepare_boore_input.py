
import numpy as np
from glob import glob
from sys import argv
from process_kikhdr import W_Cempfaultparams, hyper_epid
import os.path


#path_to_folder = argv[1]

#filenames1 = glob('/Users/jamesholt/Documents/MRes/20001006133000/20001006133000.kik/*.EW1.h')

#filenames2 = glob('/Users/jamesholt/Documents/MRes/20001006133000/20001006133000.kik/*.NS1.h')

#filenames3 = glob('/Users/jamesholt/Documents/MRes/20001006133000/20001006133000.kik/*.UD1.h')
#---------------
#SL, DDW, DDWA, RSA, SW, SWA  = W_Cempfaultparams(dip, mag)
#eqparams = np.array([h[6], h[7], h[8], srl, rw, rwa, ra, sw, swa])
#---------------
#np.savetxt('quakeparams.txt', eqparams, fmt='%10.5f', header="Quake lat(deg), Quake lon(deg), Quake Depth(km), Surface Length(km), Downdip Width(km) (with formula), Downdip W(km) (with area), Area(km^2), S width w/out area(km), S width with area(km)")

# argv in in order of argv[0] = python script name, argv[1] = full folder path, 
#... argv[2] is dip in degrees, argv[3] is the magnitude (Mw) and argv[4] 
#... is the output file name.

with open('filenames.txt', 'w') as fnames:   
    for f in filenames1:
        fnames.write(f + '\n')
    for f in filenames2:
        fnames.write(f + '\n')
    for f in filenames3:
        fnames.write(f + '\n')

def main():
    path_to_folder = argv[1]
    
    
    f1, f2, f3, f4, f5, f6 = grab_filenames(path_to_folder)
        
    write_fname_to_file(f1, f2, f3, f4, f5, f6)
    
    with open ('filenames.txt', 'r') as f:
        for fname in f:
            make_db_files(fname)


def grab_file_names(path_to_folder):
    """ This function uses glob to search for file names using the path string 
        to the folder. Glob then searches for the file extension below and 
        stores the filenames in a list. The function retuns a list for each
        of the glob searches separately. """ 
    filenames1 = glob(path_to_folder + '*.EW1.h')

    filenames2 = glob(path_to_folder + '*.NS1.h')

    filenames3 = glob(path_to_folder + '*.UD1.h')

    filenames4 = glob(path_to_folder + '*.EW2.h')
    
    filenames5 = glob(path_to_folder + '*.NS2.h')

    filenames6 = glob(path_to_folder + '*.UD2.h')

    return filenames1, filenames2, filenames3, filenames4, filenames5, filenames6

    

def write_fname_to_file(fname1, fname2, fname3, fname4, fname5, fname6):
    """ This function writes all of the file paths and names to a single file.  """
    with open('filenames.txt', 'w') as fnames:   
        for f in fname1:
            fnames.write(f + '\n')
        for f in fname2:
            fnames.write(f + '\n')
        for f in fname3:
            fnames.write(f + '\n')
        for f in fname4:
            fnames.write(f + '\n')
        for f in fname5:
            fnames.write(f + '\n')
        for f in fname6:
            fnames.write(f + '\n')

def make_db_files(fname):
    dip = float(argv[2])
    mag = float(argv[3])
 
    with open (argv[1] + argv[4], 'ab') as f:          
            
        h = np.loadtxt(fname)
        EpiD, HypD = hyper_epid(h)
        
        dbfilearray = np.array([h[10], h[11], h[12], h[23]/100, EpiD/1000, HypD/1000])
        np.savetxt(f, dbfilearray)

#def check_if_exist(path_to_file, filename):
    #if os.path.isfile(path_to_file + filename) == False:
        #break
    #else: 
        #print('Warning file already exists, please remove or rename.')
        

          
   

            
            











