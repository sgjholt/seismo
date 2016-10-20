#import sys
import numpy as np
import scipy.signal as sg
from obspy import read

#EWfile = sys.argv[2]
#NSfile = sys.argv[3]


#combine horz componants and apply calibration - *100 gives cm/s/s for strata
def combine_comps(datEW, datNS):
    """Combines the two time series """
    if min(datEW) < min(datNS):
        lowest = min(datEW)
    else:
        lowest = min(datNS)

    datEW = datEW + lowest+1
    datNS = datNS + lowest+1

    comb = (datEW**2 + datNS**2)**(1/2)

    return comb


if str(sys.argv[1]) == '-h':
   #Note - obspy gives values in m/s/s after apply calib const
   st = read(EWfile)
   st += read(NSfile)
   st.detrend()

   horz = combine_comps(st[0].data, st[1].data)*st[0].stats.calib*100
   #detrend the combined signal
   horz = sg.detrend(horz)
#Note - obspy gives values in m/s/s after apply calib const




#time = np.linspace(0,(st[0].stats.npts*st[0].stats.delta), (st[0].stats.npts))  
