# S_B.py:- script to calculate surface/borehole ratios for Kik-Net Sites
# VERSION: 
# AUTHOR(S): JAMES HOLT - UNIVERSITY OF LIVERPOOL 
# EMAIL: j.holt@liverpool.ac.uk
# 
# LAST EDIT: 
# 
# EDIT LOG: 
# 
# MAIN FEATURES:
#
# 
#
#---------------------------------modules--------------------------------------#
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from obspy.signal.konnoohmachismoothing import konno_ohmachi_smoothing
from obspy.signal.util import util_geo_km
from obspy import UTCDateTime  
import sys
import glob
import os
import gzip
from joblib import Parallel, delayed
import multiprocessing
import pandas as pd
import contextlib
from copy import deepcopy
import datetime
#---------------------------------core fns-------------------------------------#
def eventsPicker(All = True, dates=None):
    """eventsPicker will grab the directories for all events in the range 
       given by dates. By default the function will grab all directories for all
       events in the database. All is of type bool and must be set to false
       when dates are specified. dates is of type int and is inserted as a tuple
       e.g. dates = (syear,eyear,smonth,emonth). 
    """
    superDirs = []
    if All:
        [superDirs.append(
          [s.strip('.pha') for s in glob.glob(
          m+'/'+'*.pha')]) for m in Alldirs()]   

    if not All:
        syear, eyear, smonth, emonth = dates
        [superDirs.append(
          [s.strip('.pha') for s in glob.glob(
          m+'/'+'*.pha')]) for m in yearMonth(Alldirs(),syear,eyear,smonth,emonth)]

    onelist = []
    [[onelist.append(element) for element in superDirs[n]] for n in range(len(superDirs))]


    return onelist
#-----------------------------------utils--------------------------------------#

def readKiknet(fname, grabsignal=True):
    

    if fname.endswith(".gz"):
        with gzip.open(fname, 'rt', newline='\n') as f:
            if not grabsignal: # only grab headers using iter (DONT USE .readlines())
                strings = [next(f).replace("\n","").split() for x in range(15)]
            else: # grab the whole file using iter (DONT USE .readlines())
                strings = [line.replace("\r","").split() for line in f]
    else:
        with open(fname, 'rt', newline='\n') as f:
            if not grabsignal: # only grab headers using iter (DONT USE .readlines())
                strings = [next(f).replace("\n","").split() for x in range(15)]
            else: # grab the whole file using iter (DONT USE .readlines())
                strings = [line.replace("\r","").split() for line in f]

    #ASSIGN THE DATA TO APPROPRIATE VARS TO BE PASSED INTO DICT
    #station name
    stname = strings[5][2]
    #magnitude (jma)
    jmamag = float(strings[4][1])
    #frequency/dt
    freq = strings[10][2].strip('Hz')
    dt = 1/float(freq)
    #station/event lats/longs
    slat, slong, olat, olong = (float(strings[6][-1]), float(
      strings[7][-1]), float(strings[1][-1]), float(strings[2][-1]))
    #scaling factor block
    scalingF = float(strings[13][2].split('(gal)/')[0])/float(
      strings[13][2].split('(gal)/')[1]) #ARRGH WHY FORMAT LIKE THIS
    if fname.split('.')[1][-1] == "1":
        where = "Borehole: "+fname.split('.')[1][0:2]
    elif fname.split('.')[1][-1] == "2":
        where = "Surface: "+fname.split('.')[1][0:2]
    else:
        print("KNET FILE DETECTED") 
    #origin date and time (Japan)
    odate_time = UTCDateTime(datetime.datetime.strptime(strings[0][-2]+" "+strings[0][-1],'%Y/%m/%d %H:%M:%S'))-60*60*9 #-9hours for UTC time

    #recording start time (Japan)
    rdate_time = UTCDateTime(datetime.datetime.strptime(strings[9][-2]+" "+strings[9][-1],'%Y/%m/%d %H:%M:%S'))-60*60*9

    #pga
    pga = float(strings[14][-1])
    #sheight (m), eqdepth (km)
    eqdepth, sheight = (float(strings[3][-1]),float(strings[8][-1]))

    if not grabsignal: #return only the metadata
        return {"site":stname,"jmamag":jmamag,"dt":dt,"SF":scalingF,
         "origintime":odate_time,"instrument":where,"sitelatlon":(slat, slong),
         "eqlatlon":(olat, olong), "eqdepth":eqdepth, "station height":sheight,
          "pga":pga, "recordtime": rdate_time}    
    else:
        #extract the data
        data = strings[17:] #data begins at line 17 to end
        dat = np.zeros((len(data), 8)) #empty matrix of 0's to populate
        for i in range(0, len(data)):
            if len([float(l) for l in data[i]]) == 8:
        #regular expressions not needed as whitespace between numbers only
                dat[i,:] = [float(l) for l in data[i]]
            else: #append first data point (in counts) until len(array) == 8
                tmp = [float(l) for l in data[i]]
                [tmp.append(float(data[0][0])) for a in range(0, 8-len(tmp))]
                dat[i,:] = tmp

        return {"data":sig.detrend(dat.reshape(1,dat.size)[0]),"site":stname, 
         "jmamag":jmamag,"dt":dt,"SF":scalingF,"origintime":odate_time, 
         "instrument":where,"sitelatlon":(slat, slong), "eqlatlon":(olat, olong),
         "eqdepth":eqdepth, "station height":sheight,"pga":pga,
         "recordtime": rdate_time}


def calcFAS(eqDict):
    #np.fft.rfft function takes the postive side of the frequency spectrum only  
    #energy is conserved between time/freq domains
    
    #sig.tukey is a cosine taper defined in the scipy.signal package 
    # ... current taper = 5% where 0% is a square window
    FAS = np.abs(np.fft.rfft(np.pad(
      eqDict['data']*eqDict['SF'], 1000, 'constant')*sig.tukey(
      np.pad(eqDict['data']*eqDict['SF'], 1000, 'constant').size,0.05)))
    freq = np.fft.rfftfreq(len(np.pad(
      eqDict['data']*eqDict['SF'], 1000, 'constant')),eqDict['dt'])

    
    parsevalsCheck(eqDict, FAS)
    
    eqDict.update({'FAS': FAS, 'FASfreqs': freq})
    print('Added FAS/FAS-freqs to eq dictionary object.')

def parsevalsCheck(eqDict,FAS):
    t = np.sum((eqDict["data"]*eqDict["SF"])**2)
    f = np.sum(FAS**2/FAS.size)
    if abs(t-f) >= 0.01*t:
        print('Parseval Check: WARNING - energy difference is >= 1%.')
        print('Energy difference = {0}%'.format((abs(t-f)/t))*100)
    else:
        print('Parseval Check: energy is conserved.')


def Alldirs():
    pd = '/data/share/Japan/KiK-net/' #parent directory
    years = sorted(glob.glob(pd+'*')) # years in database
    months = ['01','02','03','04','05','06','07','08','09', '10', '11', '12']
    #create paths for each month for each year - make list
    alldirs = []
    for year in years:
        [alldirs.append(year+'/'+m) for m in months] 
    return alldirs

def yearMonth(directory, syear, eyear, smonth, emonth):
    subdir = []
    on = -9999
    for i in range(len(directory)):
        if int(directory[i].split('/')[5]) == syear and int(
          directory[i].split('/')[6]) == smonth:
            on = True
        if int(directory[i].split('/')[5]) == eyear and int(
          directory[i].split('/')[6]) == emonth+1:
            on = False
        if on == True:
            subdir.append(directory[i])
        if not on:
            break
    return subdir


def site_cutter(site=None, ext="EW1", Allsites = True, Dates = None, MDsearch=False, onelist=True):
    """
       site_cutter returns a list  (or list of sublists) of waveform paths 
       according to the users choice. The user can subset paths based on the 
       extension or site name. By default the function returns a list of all
       paths to all events in the database.
       
       MDsearch returns a list of directories with the extension .EW1 no matter
       what ext is set to. This will allow the user to use pythons built in set
       comparison operations to determin if data is missing. See function '
       missing_data()' for further details.
    
    """

    if site is not None:
        Allsites=False

    if not Allsites: #Check if input is correct for subset by sitename and ext
        if site is None:
            return print('Invalid Selection')
    
    if Dates is not None:
        ALL = False
    else:
        ALL = True
    Alldirs = eventsPicker(All=ALL, dates=Dates)
    
    if not Allsites and site is not None and ext is not None: #if subsetting 
        print('Cutting events by single site : {0} and ext {1}'.format(site,ext))
        
        dirs = ["{0}/{1}{2}.{3}.gz".format(
        path,site,path.split("/")[-1],ext) for path in Alldirs if os.path.isfile(
          "{0}/{1}{2}.{3}.gz".format(path,site,path.split("/")[-1],ext))]
        if not dirs:
            return print('Empty list: invalid selection.')
        else:
           return dirs
    else: #return everyfile path in every event for every site available
        Allevents = []
        if onelist:
            [[Allevents.append(event) for event in glob.glob(
              Alldirs[n]+"/*.{0}.gz".format(ext))] for n in range(len(Alldirs))];
        if not onelist: # return list with sublists separating by event
            [Allevents.append([recording for recording in glob.glob(
              Alldirs[n]+"/*.{0}.gz".format(ext))]) for n in range(len(Alldirs))];

        if MDsearch: # allows me to search if any data is missing in whole db
            return [Allevents[n].replace(ext, "EW1") for n in range(len(Allevents))];
        if not MDsearch: 
            return Allevents


def missing_data():
    """Takes the difference between the union of files with all extensions
       and subtracts the intersect to see if there are any missing data gaps for
       all available data. """
    missing_dat = list((set(site_cutter(MDsearch=True)) | set(
      site_cutter(ext="EW2",MDsearch=True)) | set(
      site_cutter(ext ="NS1",MDsearch=True)) | set(
      site_cutter(ext="NS2",MDsearch=True)) | set(
      site_cutter(ext="UD1",MDsearch=True)) | set(
      site_cutter(ext="UD2",MDsearch=True))) - (set(
      site_cutter(MDsearch=True)) & set(
      site_cutter(ext="EW2",MDsearch=True)) & set(
      site_cutter(ext ="NS1",MDsearch=True)) & set(
      site_cutter(ext="NS2",MDsearch=True)) & set(
      site_cutter(ext="UD1",MDsearch=True)) & set(
      site_cutter(ext="UD2",MDsearch=True))))


    pd = '/data/share/Japan/'
    if list(missing_dat): # if the list contains values
        with open(pd+'MissingDat.log', 'wt') as f: #write the log file
             [f.write(event+"\n") for event in list(missing_dat)]
    else:
        print('Cannot find any missing data.')

def filterMagnitude(MagRange,s=None, Ext="EW1", allsites=True):
    

    dirs = site_cutter(site=s, ext=Ext, Allsites=allsites)
    eventDicts = Parallel(n_jobs=4)(delayed(readKiknet)(event, grabsignal=False) for event in dirs)

    

    if type(MagRange) is int or type(MagRange) is float:

        if not magRangeCheck(eventDicts,MagRange)[0]:
            if MagRange > magRangeCheck(eventDicts,MagRange)[1]:
                print('Mag too big - Searching for all below')
                MagRange = magRangeCheck(eventDicts,MagRange)[1]
            if MagRange < magRangeCheck(eventDicts,MagRange)[2]:
                print('Mag too small - Searching for all equal to min mag') 
                MagRange = magRangeCheck(eventDicts,MagRange)[2]

        print("Cutting events by JMA Magnitude <= {0}".format(MagRange))
        magselect = [dirs[n] for n in range(len(
          eventDicts)) if eventDicts[n]["jmamag"] <= MagRange ]
    
    if type(MagRange) == tuple:
        if not magRangeCheck(eventDicts, MagRange)[0]:
            print('Magnitude range invalid.')
            newMagRange = [0,0]
            if MagRange[1] < magRangeCheck(eventDicts, MagRange)[2]:
                newMagRange[1] = magRangeCheck(eventDicts, MagRange)[2]
            else: 
                newMagRange[1] = MagRange[1]
            if MagRange[0] > magRangeCheck(eventDicts, MagRange)[1]:
                newMagRange[0] = magRangeCheck(eventDicts, MagRange)[1]
            else: 
                newMagRange[0] = MagRange[0]
            del MagRange
            MagRange = tuple(newMagRange)

        print("Cutting events by JMA Magnitude {1} <= Mag <= {0}".format(
          MagRange[0], MagRange[1]))
        lowerthan = [dirs[n] for n in range(
          len(eventDicts)) if eventDicts[n]["jmamag"] <= MagRange[0]]
        greaterthan = [dirs[n] for n in range(
          len(eventDicts)) if eventDicts[n]["jmamag"] >= MagRange[1]]
        magselect = list(set(lowerthan) & set(greaterthan))
        
    if not magselect:
        print('Magnitude range not available')

    return magselect

def event_tables(oneCORE = False, NCORE=None):
    allEvents = site_cutter(onelist=False)
    maxCPUs = multiprocessing.cpu_count()
    if NCORE > maxCPUs or NCORE is None:
        NCORE = maxCPUs
        print('Not enough cores, setting n_jobs to {0}'.format(maxCPUs))
    Parallel(n_jobs=int(NCORE))(delayed(useful_metadata)(
      eventList) for eventList in allEvents)
    if not oneCORE:
        [useful_metadata(eventList) for eventList in allEvents]; #
  

def useful_metadata(eventDirs, save=True):
    eventDirs_copy = deepcopy(eventDirs) #dont want to modify input list

    wfDict = [readKiknet(wf, grabsignal=False) for wf in eventDirs];
    wfDict += [readKiknet(wf.replace(
      "EW1", "EW2"),grabsignal=False) for wf in eventDirs]
    eventDirs_copy *= 2
    eventDirs_copy = [Dir.replace(".EW1.gz", "") for Dir in eventDirs_copy]
    table = pd.DataFrame(columns = ("site", "jmamag", "pga", "site_lat", "site_lon",
      "station_height", "instrument", "eq_lat", "eq_lon", "eq_depth","Repi", 
      "Rhypo", "o_year", "o_month", "o_day", "o_time","rec_start_time", "eventid","path" ))
    n = 1
    for wf in wfDict:
        table.loc[n] = wf["site"], wf["jmamag"], wf["pga"],wf["sitelatlon"][0], wf[
          "sitelatlon"][1], wf["station height"],wf["instrument"].split(
          ':')[0], wf["eqlatlon"][0],wf["eqlatlon"][1], wf["eqdepth"],calcEpiHypo(
          wf)[0],calcEpiHypo(wf)[1], wf["origintime"].year, wf[
          "origintime"].month, wf["origintime"].day, wf["origintime"].strftime(
          '%H:%M:%S'), wf["recordtime"], eventDirs_copy[0].split(
          '/')[-2],eventDirs_copy[n-1]
        n += 1

    if save: #should save by default as likely this is what you want
        silent_remove(eventDirs[0].replace(eventDirs[0].split(
          '/')[-1], "")+eventDirs[0].split('/')[-2]+"_db.csv")
            
        table.to_csv(eventDirs[0].replace(eventDirs[0].split(
          '/')[-1], "")+eventDirs[0].split('/')[-2]+"_db.csv")
        
        print("Saved database as db_{0}.csv for event {0}".format(
          eventDirs[0].split('/')[-2]))
       
            
    else: # save = False, therefore return table as pd.DataFrame object.
        return table
    
    


def calcEpiHypo(wf):
    """
       Uses Obspy's distance converter tool to calulate the distance delta
       between site and eq; then converts to km from deci degs. 
       Takes in waveform dictionary object returned by the readKiknet() 
       function. Better to use Obspy's func because it uses elliptical earth 
       model = accurate distances.
    """
    #extract relevent params
    slat, slon = wf["sitelatlon"]
    eqlat, eqlon = wf["eqlatlon"]
    eqdepth = wf["eqdepth"] 
    sheight = wf["station height"] / 1000 #convert to km
    #do the delta and conversion to km
    dx,dy = util_geo_km(eqlon,eqlat,slon,slat)
    #calc Repi (dx**2 + dy**2)**0.5 and Rhyp (dx**2 + dy**2 + dz**2)**0.5  
    Repi = np.sqrt(dx**2 + dy**2)
    Rhyp = np.sqrt(dx**2 + dy**2 + (eqdepth+sheight)**2)

    return (Repi,Rhyp)

       
def magRangeCheck(eventDicts, MagRange):
    MaxMag = max([eventDicts[n]["jmamag"] for n in range(len(eventDicts))])
    MinMag = min([eventDicts[n]["jmamag"] for n in range(len(eventDicts))])
    if type(MagRange) is float or type(MagRange) is int:
        if MagRange < MinMag or MagRange > MaxMag:
            return False, MaxMag, MinMag
        else:
            return True, MaxMag, MinMag
    elif type(MagRange) is tuple:
        if MagRange[0] > MaxMag or MagRange[1] < MinMag:
            return False, MaxMag, MinMag
        else:
            return True, MaxMag, MinMag

def silent_remove(filename):
    """
       The silent_remove() function uses the os.remove function to remove a 
       given file. The contextlib module suppresses the error given if the file 
       does not exist, allowing the script to continue executing. 


       USAGE: silent_remove(path+filename)
       parent_dir = /home/directory_of_interest/
       fname = "file.txt"
       E.G. - silent_remove(parent_dir + fname) 
    """
    with contextlib.suppress(FileNotFoundError):
        os.remove(filename)
#---------------------------------run cmds-------------------------------------#
if __name__ = '__main__'
    event_tables(NCORE=4)

