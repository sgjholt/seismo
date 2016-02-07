import numpy as np
import os
from glob import glob
from build_GMdb import grab_file_names

def chkEmptyhdr(fname):
    for filename in fname:
        f = np.loadtxt(filename)

        if f.size > 0:
            continue
        

        else:
            print('File {} is empty, please check the source file.'.format(filename))
    print('Check is complete.')   

  
    
