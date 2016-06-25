#OBSPY UTILITIES


#-------------------------------------------------------------------------------
#IMPORT MODULES
import numpy as np
from glob import glob

#-------------------------------------------------------------------------------


def GetFiles(surf_bore_all, path):
    """ THIS FUNCTION RETURNS THE FILE LIST SPECIFIED BY THE USER.
        USAGE: GetFiles('surface', '20001006133000.kik') - returns a python list
        object containing the file names relating the to the surface seismometers
        of the kik-net network. """

    path = '/Volumes/J_Holt_HDD/MRes/Modules/Thesis/Data/' + str(path) + '/'

    if surf_bore_all == 'surface':
        flist = glob(path + '*.*2')

    elif surf_bore_all == 'borehole':
        flist = glob(path + '*.*1')

    else:
        flist = glob(path + '*.*')
        flist = [f.replace('.h', '') for f in flist]

    return flist


 

