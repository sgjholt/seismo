from obspy.core import read
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from utils import grab_file_names
import scipy.signal as sg
import glob



#-------------------------------------------------------------------------------
#argv tracker
#argv[0] - script name
#argv[1] - eq folder path (eg. '/Users/jamesholt/Desktop/20160416012500.kik')
#argv[2] - eq database file name (eg. '20160415_M7_DB_complete.txt')
#argv[3] - surface or downhole? (eg. 'Surface' / 'Downhole')
#-------------------------------------------------------------------------------

#freq = dat[:,0]
#FAS = dat[:,1]

#plt.title(
#'Fourier Acceleration Spectrum for Mw 6.6 Earthquake\n Station:GNMH08 - EW component',fontsize=16, y=1.01) 
#plt.loglog(freq[(freq>0.04) & (freq<=10)], FAS[(freq>0.04)&(freq<=10)], 'b-')
#plt.xlabel('FREQUENCY [HZ]',fontsize=14)
#plt.ylabel('FAS [M/S]',fontsize=14)
#plt.savefig(
#'/home/james/Dropbox/MRes/Modules/Thesis/poster/FAS.pdf')
#plt.close()






def main(argv):
    
    st = streamBuild(argv[1], argv[2]) #load in data stream
    s_pick = S_picker(st) #pick the S-wave arrivals
    s_window(st, s_pick) #calculate the positions of window
    print("Grabbing S waveforms.")
    WindowGrabber(st) #grab the s-windows in time domain
    print("Generating FAS.")
    FASmaker(st, argv)
    num = [0.3, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    #for i in range(0, int(len(num))):
        
    #    freqPicker(st, num[i], argv)
    #argv[3] = 'Surface'
    for i in range(0, int(len(num))):
        
        argv = freqPicker(st, num[i], argv)

    return st



def Freqy(List, argv):
    st = streamBuild(argv[1], argv[2])
    for i in range(0, int(len(List))):
        FAS = np.loadtxt(List[i])
        st[i].stats["FAS"] = {}
        st[i].stats["FAS"]["Spectrum"] = FAS[:,1]
        st[i].stats["FAS"]["Frequency"] = FAS[:,0]
       ls
    return st


def listmaker():
    List = glob.glob('*.EW1.FAS')
    List += glob.glob('*.NS1.FAS')
    List += glob.glob('*.UD1.FAS')
    List += glob.glob('*.EW2.FAS')
    List += glob.glob('*.NS2.FAS')
    List += glob.glob('*.UD2.FAS')
    return List

def multiQuakePlot(List, FREQ, argv, on_off):
    if argv[3] == 'Downhole':
        F = 'FAS_dist_DWN_'+str(FREQ)+'.txt'
    else:
        F = 'FAS_dist_SRF_'+str(FREQ)+'.txt'
    for f in List:
        DAT = np.loadtxt('/Volumes/J_Holt_HDD/MRes/Modules/Thesis/Data/'+f+'/'+F)
        #print(len(DAT[1]))
        Mo = np.loadtxt(
        '/Volumes/J_Holt_HDD/MRes/Modules/Thesis/Data/'+f+'/'+'Mo.txt')
        SF = FASscaler(FREQ, 7, Mo)
        ax1 = plt.subplot(121)
        plt.ylabel('FAS [M/S]')
        plt.xlabel('RJB [KM]')
        plt.loglog(DAT[1]/1000, DAT[0],'.')
        plt.setp(ax1.get_xticklabels(), fontsize=14)
        plt.setp(ax1.get_yticklabels(), fontsize=14)
        ax1.set_axis_bgcolor('w')
        if argv[3] == 'Downhole':
            plt.suptitle(
            'FAS @ '+str(
             FREQ)+' HZ Raw [left] and Shifted to Mw 7 [right] : Borehole', fontsize=16)
        else:
            plt.suptitle(
            'FAS @ '+str(
            FREQ)+' HZ Raw [left] and Shifted to Mw 7 [right] : Surface', fontsize=16)
        ax2 = plt.subplot(122, sharex = ax1, sharey=ax1)
        plt.loglog(DAT[1]/1000, DAT[0]*SF,'.')
        plt.xlabel('RJB [KM]')
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), fontsize=14)
        ax2.set_axis_bgcolor('w')
    if on_off == 'on':
        #smooth_plotter(List, FREQ, 30, argv)
        plt.show()
    else:
        if argv[3] == 'Downhole':
            plt.savefig(
            '/Users/jamesholt/Dropbox/MRes/Modules/Thesis/poster/FAS_DWN_'+str(FREQ)+'_HZ.pdf')
            plt.close() #DONT FORGET TO CLOSE THE PLOT IF SAVING FIGURES IN BATCH!!num
        else:
            plt.savefig(
            '/Users/jamesholt/Dropbox/MRes/Modules/Thesis/poster/FAS_SRF_'+str(FREQ)+'_HZ.pdf')
            plt.close()



def smooth_plotter(List, FREQ, binNo, argv, on_off):
    #collect all the data into one vector for scaled and not scaled data
    out, out_SF = collectDATA(List, FREQ, argv) 
    #bin the data appropriately, include information about the data variance
    stats = binner(out, binNo)
    statsSF = binner(out_SF, binNo)
    #calculate the difference between the max/min and median points for variance plot
    #(errorbars)
    dat_range = np.array([(stats[5]), (stats[5])])
    dat_rangeSF = np.array([(statsSF[5]), (statsSF[5])])
    #plot the log of FAS vs distance for binned data
    plt.plot(
    stats[0]/1000, stats[3], 'rs', markersize=5, alpha = 0.7, label = 'Raw Data' )
    plt.plot(
    stats[0]/1000, statsSF[3], 'bs', markersize=5, label = 'Shifted to Mw 7')
    red_sq = mlines.Line2D([], [], color='red', marker='s',
                          markersize=15, label='Raw Data')   
    blue_sq = mlines.Line2D([], [], color='blue', marker='s',
                          markersize=15, label='Shifted to Mw 7')

    plt.legend(handles=[red_sq, blue_sq], loc='lower left')

    plt.errorbar(stats[0]/1000, statsSF[3], dat_rangeSF, fmt = 'bs', capsize=10)
    plt.errorbar(stats[0]/1000, stats[3], dat_range, fmt = 'rs', capsize=0, linewidth=5, alpha=0.3)
    #errorfill(stats[0]/1000, stats[3], dat_range, color = 'r')
    #errorfill(stats[0]/1000, statsSF[3], dat_rangeSF, color = 'b')
    #labels and stuff
    plt.xlabel('RJB [KM]', fontsize=14)
    plt.ylabel('LOG10(FAS) [M/S]', fontsize=14)
    
    if argv[3] == 'Downhole':
        plt.title(
        'Smoothed Plot of the Median of \n'+str(
         binNo)+' Data Bins for FAS @ '+str(
         FREQ)+' HZ : Borehole', y=1.01, fontsize=16) 
    else:
        plt.title(
        'Smoothed Plot of the Median of \n'+str(
         binNo)+' Data Bins for FAS @ '+str(
         FREQ)+' HZ : Surface', y=1.01, fontsize=16)
    if on_off == 'on':
        #smooth_plotter(List, FREQ, 30, argv)
        plt.show()
    else:
        if argv[3] == 'Downhole':
            plt.savefig(
            '/Users/jamesholt/Dropbox/MRes/Modules/Thesis/poster/smooth_FAS_DWN_'+str(
            FREQ)+'_HZ_bins='+str(binNo)+'.pdf')
            plt.close()
            
        else:
            plt.savefig(
            '//Users/jamesholt/Dropbox/MRes/Modules/Thesis/poster/smooth_FAS_SRF_'+str(
             FREQ)+'_HZ_bins='+str(binNo)+'.pdf')
            plt.close()

def errorfill(x, y, yerr, color=None, alpha_fill=0.3, ax=None):
    ax = ax if ax is not None else plt.gca()
    if color is None:
        color = ax._get_lines.color_cycle.next()
    if np.isscalar(yerr) or len(yerr) == len(y):
        ymin = y - yerr
        ymax = y + yerr
    elif len(yerr) == 2:
        ymin = yerr[0] 
        ymax = yerr[1]
    ax.plot(x, y, color=color)
    ax.fill_between(x, ymax, ymin, color=color, alpha=alpha_fill)



def binner(dat, noBins):
    dist = dat[1]  #distance data
    data = dat[0]  #actual data
    hist = np.histogram(dist, noBins)
    bins = hist[1] 
    stats  = np.zeros([6, (int(len(bins))-1)])
    
    
    for i in range(1, int(len(bins))):
        if i == 1:
            bin1 = np.array([ data[dist < bins[i]], dist[dist < bins[i]] ])
        if i > 1: 
            bin1 = np.array(
            [data[(dist >=bins[i-1]
            ) & (dist<bins[i])], dist[(dist >=bins[i-1]) & (dist<bins[i])]])
        #meandist = bin1[1].mean()
        #binmaxdist = bin1[1].max()
        #binmindist = bin1[1].min()
        #meanpt = bin1[0].mean()
        #binmax = bin1[0].max()
        #binmin = bin1[0].min()
        #stats[:,i-1] = np.array(
        #[np.median(bin1[1]), bin1[1].max(),bin1[1].min(),np.mean(bin1[0]),bin1[0].max(),bin1[0].min()]) 


        meanbin = np.mean(np.log10(bin1[0]))
        varbin = np.sum(((np.log10(bin1[0])-meanbin)**2))/len(bin1[0])
        std = varbin**(1/2)
        stats[:,i-1] = np.array(
        [np.median(bin1[1]), bin1[1].max(),bin1[1].min(),meanbin,varbin,std]) 
    return stats 
        
            
   #bin2 = np.array([data[(dist >=50) & (dist<100)], dist[(dist >=50) & (dist<100)]]) 

   
    


def collectDATA(List, FREQ, argv):
    if argv[3] == 'Downhole':
        F = 'FAS_dist_DWN_'+str(FREQ)+'.txt'
    else:
        F = 'FAS_dist_SRF_'+str(FREQ)+'.txt'
   
    for i in range(0, int(len(List))):
        #if i > 1:
            #OUTVEC_SF = np.concatenate((OUTVEC_SF,tmp_SF), axis=1)
            #OUTVEC = np.concatenate((OUTVEC,tmp), axis=1)
        f = str(List[i])
        DAT = np.loadtxt('/Volumes/J_Holt_HDD/MRes/Modules/Thesis/Data/'+f+'/'+F)

        Mo = np.loadtxt(
        '/Volumes/J_Holt_HDD/MRes/Modules/Thesis/Data/'+f+'/'+'Mo.txt')
        SF = FASscaler(FREQ, 7, Mo)
        tmp_SF = np.array([(DAT[0]*SF), (DAT[1])])

        tmp = DAT

        if i == 0:
            OUTVEC_SF = tmp_SF
            OUTVEC = tmp
            #print('Beginning with {0} length {1}'.format(f, len(tmp[1])))
        if i > 0:
            OUTVEC_SF = np.concatenate((OUTVEC_SF,tmp_SF), axis=1)
            OUTVEC = np.concatenate((OUTVEC,tmp), axis=1)
            #print('Concatenating to vector {0} length {1}'.format(f, len(tmp[1])))
        #print('final length = {0}'.format(len(OUTVEC[1])))
        
    return OUTVEC, OUTVEC_SF
        
def histogram(List, FREQ, argv, binwidth):
    dat, datSF = collectDATA(List, FREQ, argv)
    plt.figure(figsize=(14, 8))
    plt.suptitle('Histogram of FAS before Brune source correction \n [left] and after [right] for ' +str(FREQ)+ 'Hz :' + str(argv[3]), fontsize=16, y=0.99)
    
    plt.subplot(121)
    plt.xlabel('FAS')
    plt.ylabel('No of Occurrences')
    
    plt.hist(np.hstack(dat[0]), bins=np.arange(min(dat[0]), max(dat[0]) + binwidth, binwidth))
    plt.subplot(122)
    plt.xlabel('FAS')
    
    plt.hist(np.hstack(datSF[0]), bins=np.arange(min(datSF[0]), max(datSF[0]) + binwidth, binwidth))
   
    plt.show()    
        



def plotter(List, num):

    i = 0
    for f in List:
    
        ratio = FASscaler(num[i], 7, 4.42e19)
        DAT = np.loadtxt(f)
        plt.loglog(DAT[1]/1000, DAT[0]*ratio, '.')
        i += 1
    plt.show()

def FASscaler(f, Mw, MoREAL):
    """Calculates the scaling factor to shift FAS to that of a reference 
    earthquake. This is the ratio between a reference Brune source at reference 
    magnitude and the real event. Requres input of the real Seismic Moment, 
    and frequency.  """   
    Mo = idealMo(Mw)
    C = (0.4906 * 3500)
    fcReal = C*(10E6/MoREAL)**(1/3)
    fc = C*(10E6/Mo)**(1/3)
    ratio = (Mo*(1+(f/fcReal)**2))/(MoREAL * (1+(f/fc)**2))
    return ratio


def idealMo(momentMAG): 
    """Using Kanamori's formula, calculates the ideal Mo for a given Mw. """ 
    Mo = 10**( ( momentMAG*1.5 ) + (6.03*1.5) )
    return Mo

def freqPicker(st, FREQ, argv):

    halfst = int(len(st)/2)
    group = int(len(st)/6)
    value = np.zeros(group)
    dist = np.zeros(group)
    for i in range(0, int(len(st)/6)):
    
        if argv[3] == 'Downhole':
            FASew = st[i].stats.FAS.Spectrum
            FASns = st[i+group].stats.FAS.Spectrum
            FASud = st[i+group*2].stats.FAS.Spectrum
            freq = st[i].stats.FAS.Frequency
            d = st[i].stats.distance
        else:
            FASew = st[i+halfst].stats.FAS.Spectrum
            FASns = st[i+halfst+group].stats.FAS.Spectrum
            FASud = st[i+halfst+group*2].stats.FAS.Spectrum
            freq = st[i+halfst].stats.FAS.Frequency
            d = st[i+halfst].stats.distance

        triFAS = geoMean(FASew, FASns, FASud)
        nearest = find_nearest(freq, FREQ)    
        where = np.where(freq == nearest)

        value[i] = triFAS[where]       
        dist[i] = d
    if argv[3] == 'Downhole':
        np.savetxt((argv[1]+"/FAS_dist_DWN_"+str(FREQ)+".txt"), (value, dist))
    else: 
        np.savetxt((argv[1]+"/FAS_dist_SRF_"+str(FREQ)+".txt"), (value, dist))

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def FASmaker(st, argv):
    
    for i in range(0, len(st)):
        FAS = calcFAS(st, st[i].stats.swindow, i)        
        st[i].stats["FAS"] = {}
        st[i].stats["FAS"]["Spectrum"] = FAS[:,1]
        st[i].stats["FAS"]["Frequency"] = FAS[:,0]
        np.savetxt(argv[1]+"S_WindFAS_"+str(st[i].stats.station)+"."+str(
        st[i].stats.channel)+".FAS", FAS)
        print('FAS {0} out of {1}'.format(i+1, int(len(st)))) 


def geoMean(file1, file2, file3):
    """Calculates the geometric mean of three spectral componants (e.g. N-S and
    E-W, U-D componants. The files must already be loaded into python using np.load.
    USAGE: geomean = geoMean(EW, NS). """

    geoMean = (file1*file2*file3)**(1/3)

    return geoMean



def calcFAS(st, i):
    
    n = len(st[i].data) #length of the signal
    freq = np.fft.rfftfreq(n, st[i].stats.delta) #calculate frequencies 
    wind = np.blackman(n) #blackman window to reduce spectral leaks
    #calculate fas - rfft because real signals have hermitian symmetry.
    fas = np.abs(np.fft.rfft(st[i].data*st[i].stats.calib*wind)) / (n/2) 
    FAS = [freq[3:int(len(freq)-1)], fas[3:int(len(freq)-1)]] 
    return np.transpose(FAS)

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




def S_picker(st):
    """S_picker requires an input of an obspy stream variable to function.
        This is a simple algorithm which seeks the maximum amplitude of each 
        waveform (for NS/EW - Surface/Borehole) and takes the average position 
        as the S-wave pick. Returns the pick in seconds. """
    S_pick = np.zeros(int(len(st)/6))
    for i in range(0, int(len(st)/6)):
        Spt_EW_bore = np.where(st[i].data == np.max(st[i].data)) #EW1 borehole
        
        Spt_NS_bore = np.where(st[i+int(len(st)/6)].data == np.max(
        st[i+int(len(st)/6)].data)) #NS1 borehole

        Spt_EW_surf = np.where(st[i+int(len(st)/2)].data == np.max(
        st[i+int(len(st)/2)].data)) #EW2 surface
        
        Spt_NS_surf = np.where(st[i+int(len(st)/6)+int(len(st)/2)].data == np.max(
        st[i+int(len(st)/6)+int(len(st)/2)].data)) #NS2 surface   

        meanSpt = (
        Spt_EW_bore[0] + Spt_NS_bore[0] + Spt_EW_surf[0] + Spt_NS_surf[0])/4
       
            
        pick = meanSpt*st[i].stats.delta
        
        S_pick[i] = np.floor(pick[0])        
    return S_pick


def s_window(st, s_pick):
    #Max = maxNPTS(st)
    #s_winds = np.zeros([int(len(s_pick)), int(Max)+1])
    halfG = int(len(st)/2)
    compG = int(len(st)/6)
    for i in range(0, int(len(s_pick))):

        if s_pick[i] <= 1.0:
            s_window = np.arange(
            s_pick[i], s_pick[i]+(
            st[i].stats.npts*st[i].stats.delta-s_pick[i])+(
            st[i].stats.delta)-1, st[i].stats.delta)
        elif s_pick[i] < 5.0 and s_pick[i] > 1.0:
            s_window = np.arange(
            s_pick[i]-1, s_pick[i]+(
            st[i].stats.npts*st[i].stats.delta-s_pick[i])+(
            st[i].stats.delta)-1, st[i].stats.delta)
        else:
            s_window = np.arange(
            s_pick[i]-5, s_pick[i]+(
            st[i].stats.npts*st[i].stats.delta-s_pick[i])+(
            st[i].stats.delta)-1, st[i].stats.delta)      
            
        st[i].stats["swindow"] = s_window/st[0].stats.delta
        st[i+compG].stats["swindow"] = s_window/st[0].stats.delta
        st[i+compG*2].stats["swindow"] = s_window/st[0].stats.delta
        st[i+halfG].stats["swindow"] = s_window/st[0].stats.delta
        st[i+compG+halfG].stats["swindow"] = s_window/st[0].stats.delta
        st[i+compG*2+halfG].stats["swindow"] = s_window/st[0].stats.delta
        #Len = Max - len(s_window)
        #s_window = np.lib.pad(
        #s_window, (0, int(Len)+1), 'constant', constant_values=(0))
        
        #s_winds[i,:] = s_window[0:(int(Max)+1)] 

    #return s_winds/st[0].stats.delta


def WindowGrabber(st):
    #s_winds = getBIG(s_winds) #make it bigger to grab all componants
   
    for n in range(0, int(len(st))):   
        dat = st[n].data
        s_wind = st[n].stats.swindow
        try:        
            for i in range(0, len(s_wind)):
                s_wind[i] = dat[s_wind[i]]

        except: IndexError 
               
        st[n].stats["swindow"] = s_wind*st[n].stats.calib #for real units
   
    

def getBIG(mat):
    a = mat
    b = a
    a = np.concatenate((a, b), axis=0)
    a = np.concatenate((a, b), axis=0)
    
    return a

def maxNPTS(st):
    npts = np.zeros(int(len(st)/6))
    for i in range(0, int(len(st)/6)):
        npts[i] = st[i].stats.npts
        Max = max(npts)
    return Max

#if __name__ == __main__:
    #main()
 
