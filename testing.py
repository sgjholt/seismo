from obspy.core import read
import numpy as np
import matplotlib.pyplot as plt
from obspy.signal.trigger import ar_pick
from utils import grab_file_names


path = '/Users/jamesholt/Desktop/20160416012500.kik/'

def PS_AR_Picker(stream, flag):
    
    group_size = int(len(stream)/6) #amount of stations in each group ...
                                    #(EW1, EW2 etc..)
    half_group = int(len(stream)/2)

    if flag == 'Downhole':     
        for i in range(0, group_size):
            tr1 = st[i].filter(
                'lowpass', freq=1.0, corners=2, zerophase=True)
            tr2 = st[i+group_size].filter(
                'lowpass', freq=1.0, corners=2, zerophase=True)
            tr3 = st[i+group_size*2].filter(
                'lowpass', freq=1.0, corners=2, zerophase=True)
            df = tr1.stats.sampling_rate
            p_pick, s_pick = ar_pick(
            tr1.data, tr2.data, tr3.data, df, 1.0, 20.0, 1.0, 0.1, 4.0, 1.0, 2, 8, 0.1, 0.2)
            print('{}s {}s from station {}' .format(p_pick, s_pick, tr1.stats.station))
    if flag == 'Surface':
        for i in range(0, group_size):
            tr1 = st[i].filter(
                'lowpass', freq=1.0, corners=2, zerophase=True)
            tr2 = st[half_group+i+group_size].filter(
                'lowpass', freq=1.0, corners=2, zerophase=True)   
            tr3 = st[half_group+i+group_size*2].filter(
                'lowpass', freq=1.0, corners=2, zerophase=True)
            df = tr1.stats.sampling_rate
            p_pick, s_pick = ar_pick(
            tr1.data, tr2.data, tr3.data, df, 1.0, 20.0, 1.0, 0.1, 4.0, 1.0, 2, 8, 0.1, 0.2)
            print('{}  s {}  s from {}' .format(p_pick, s_pick, tr1.stats.station))

st[0].stats["coordinates"] = {} # add the coordinates to your dictionary,
                                  #needed for the section plot
st[0].stats["coordinates"]["latitude"] = st[0].stats.knet.stla
st[0].stats["coordinates"]["longitude"] = st[0].stats.knet.stlo
st[0].stats["elevation"] = st[0].stats.knet.stel

for i in range(0, len(st)):
    st[i].stats["coordinates"] = {} # add the coordinates to your dictionary,
    st[i].stats["coordinates"]["latitude"] = st[i].stats.knet.stla
    st[i].stats["coordinates"]["longitude"] = st[i].stats.knet.stlo
    st[i].stats["elevation"] = st[i].stats.knet.stel

df = tr1.stats.sampling_rate
p_pick, s_pick = ar_pick(tr1.data, tr2.data, tr3.data, df, 1.0, 20.0, 1.0, 0.1, 4.0, 1.0, 2, 8, 0.1, 0.2)

def streamBuild(path):
    fnames = grab_file_names(path, 2)
    for i in range(0, len(fnames)):
        if i == 0:
            st = read(fnames[i])
        if i > 0:
            st += read(fnames[i])

    return st












