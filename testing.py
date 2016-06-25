from obspy.core import read
import numpy as np
import matplotlib.pyplot as plt
from obspy.signal.trigger import ar_pick
from utils import grab_file_names
from obspy.taup import TauPyModel

path = '/Users/jamesholt/Desktop/20160416012500.kik/'
db = '20160415_M7_DB_complete.txt'

def PS_AR_Picker(st, flag, lower_F, upper_F):
    upper_F = float(upper_F)
    
    group_size = int(len(st)/6) #amount of stations in each group ...
                                    #(EW1, EW2 etc..)
    half_group = int(len(st)/2)

    
    st.taper(max_percentage=0.1)
    st.filter(
    'bandpass', freqmin = lower_F, freqmax=upper_F, corners=2, zerophase=True)

    if flag == 'Downhole':     
        for i in range(0, group_size):
            tr1 = st[i]
            tr2 = st[i+group_size]
            tr3 = st[i+group_size*2]

            df = tr1.stats.sampling_rate
            p_pick, s_pick = ar_pick(
            tr1.data, tr2.data, tr3.data, df, 1.0, 20.0, 1.0, 0.1, 4.0, 1.0, 2, 8, 0.1, 0.2)
            print('{}s {}s from station {}' .format(p_pick, s_pick, i+1))
    if flag == 'Surface':
        for i in range(0, group_size):
            tr1 = st[i]
            tr2 = st[half_group+i+group_size] 
            tr3 = st[half_group+i+group_size*2]

            df = tr1.stats.sampling_rate
            p_pick, s_pick = ar_pick(
            tr1.data, tr2.data, tr3.data, df, 1.0, 20.0, 1.0, 0.1, 4.0, 1.0, 2, 8, 0.1, 0.2)
            print('{}  s {}  s from {}' .format(p_pick, s_pick, i+half_group+1))



def streamBuild(path,db):
    fnames = grab_file_names(path, 2)
    for i in range(0, len(fnames)):
        if i == 0:
            st = read(fnames[i])
        if i > 0:
            st += read(fnames[i])
    for i in range(0, len(st)):
        st[i].stats["coordinates"] = {} # add the coordinates to your dict
        st[i].stats["coordinates"]["latitude"] = st[i].stats.knet.stla
        st[i].stats["coordinates"]["longitude"] = st[i].stats.knet.stlo
        st[i].stats["elevation"] = st[i].stats.knet.stel
    
    datB = np.loadtxt(path+db, skiprows=1)
    datB = datB[:,8] #rjb
    
    for i in range(0, len(datB)):
        st[i].stats["distance"] = datB[i]*1000 #dist in meters
    st.detrend()
    
    return st



def tauModelCheck(st, i)    
epiD = ((st[i].stats.knet.evla - st[i].stats.coordinates.latitude)**2 + (st[i].stats.knet.evlo - st[i].stats.coordinates.longitude)**2)**0.5 
    
model = TauPyModel(model="iasp91")
arrivals = model.get_travel_times(source_depth_in_km=st[i].stats.knet.evdp, distance_in_degree=epiD, phase_list=["P", "S", "p", "s"])
    
   


