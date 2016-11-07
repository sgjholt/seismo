# conv4strata - A SCRIPT TO CONVERT KIK-NET DATA FOR USE IN THE STRATA PACKAGE
# VERSION: 1.1.1
# AUTHOR: JAMES HOLT - UNIVERSITY OF LIVERPOOL 
# EMAIL: j.holt@liverpool.ac.uk
# LAST EDIT: 2016/11/07
# 
# 2016/11/07: 
# Added: Ability to parse gzipped files as well as regular ASCII files.
# Improved algorithm effeciency when parsing the read list (less steps).
# MAIN FEATURES
#
# DATA PREVIEW:

# 1. TAKES ANY SINGLE COMPONANT KIK-NET FILE (E.G. 'FOO.EW1') AND DISPLAYS 
#  DISPLAYS THE CORRESPONDING WAVEFORM  

# 2.TAKES ANY TWO COMPONANT KIK-NET FILES (E.G. 'FOO.EW1') AND WILL COMBINE THE 
#  TWO RECORDS (AS MAGNITUDE OF VECTORS IN TIME DOMAIN) AND DISPLAYS THE 
#  CORRESPONDING WAVEFORM.  
 
# FILE CONVERSION:

# 1.TAKES ANY SINGLE COMPONANT KIK-NET FILE (E.G. 'FOO.EW1') AND WILL CONVERT 
#  THIS DATA TO PEER.AT2 FORMAT (STYLE).

# 2.TAKES ANY TWO COMPONANT KIK-NET FILES (E.G. 'FOO.EW1') AND WILL COMBINE THE 
#  TWO RECORDS (AS MAGNITUDE OF VECTORS IN TIME DOMAIN) THEN CONVERT THE DATA TO 
#  A SINGLE PEER.AT2 FORMAT (STYLE) FILE.

# PLEASE NOTE PEER.AT2 FORMATS USUALLY HAVE ONLY 5 COLUMNS - THIS DOES NOT 
# CONVERT THE DATA INTO A (Npts/5,5) STRUCTURE BUT RETAINS THE (Npts/8, 8) 
# STRUCTURE AS THIS IS ADEQUATE FOR USE IN STRATA. 


#---------------------------------modules--------------------------------------#
import sys
import numpy as np
import scipy.signal as sg
#from obspy import read
import matplotlib.pyplot as plt
#import re
import os
import gzip

#---------------------------------core fns-------------------------------------#
def preview_plot():

    fname = str(sys.argv[2]) 
       
    dat,stname,mag,dt,sf,odate = read_kiknet(fname) 
   
    npts = len(dat)*8
    
    time = np.arange(0, npts*dt,dt) #subtract dt to get len(npts) for time
    

    plt.plot(time, sg.detrend(dat.reshape(1, npts)[0])*sf, 'k-')
    plt.title('Mag: '+str(mag)+' JMA.'+' Station: '+stname+(
        ' Origin Time/Date: ') +odate+' (Japan)')
    plt.xlabel('time[s]')
    plt.ylabel('acc [cm/s/s]')
    plt.show() 
    print('Program Terminated')


def preview_horz():
    
    fname1 = str(sys.argv[2])
    fname2 = str(sys.argv[3])

    data1,stname,mag,dt,sf,odate = read_kiknet(fname1)
    data2,stname,mag,dt,sf,odate = read_kiknet(fname2)
    
    horz = combine_comps(data1, data2)

    npts = len(horz)*8
    
    time = np.arange(0, npts*dt,dt) #subtract dt to get len(npts) for time
   
    plt.plot(time, horz.reshape(1, len(horz)*8)[0]*sf, 'k-')
    plt.title('Combined Horizontal. Mag: '+str(mag)+' JMA.'+' Station: '+stname+(
        ' Origin Time/Date: ') +odate+' (Japan)')
    plt.xlabel('time[s]')
    plt.ylabel('acc [cm/s/s]')
    plt.show() 
    print('Program Terminated')

def simple_conv():
    
    fname = str(sys.argv[2]) 
    dat,stname,mag,dt,sf,odate = read_kiknet(fname) 
    
    npts = len(dat)*8

      
    convert2Strata(fname, sg.detrend(dat.reshape(1, npts)).reshape(len(dat),8)*sf, dt, npts)

def horz_conv():

    fname1 = str(sys.argv[2])
    fname2 = str(sys.argv[3])

    data1,stname,mag,dt,sf,odate = read_kiknet(fname1)
    data2,stname,mag,dt,sf,odate = read_kiknet(fname2)
    
    horz = combine_comps(data1, data2)

    npts = len(horz)*8  

    fname1 = fname1.split('.')
    fname = fname1[0]+'.horz' 

    convert2Strata(fname, horz*sf, dt, npts) 


#-----------------------------------utils--------------------------------------#
#parse the kiknet file - return python friendly outputs
def convert2Strata(fname, dat, dt, npts):
    """Takes array in numpy array format and saves text file in the PEER .AT2
    format for Strata. An example file of the PEER format can be found in the
    'examples' folder for Strata. Data must already be detrended and have 
    factor applied."""
    fname = fname.split('.')

    np.set_printoptions(formatter={'float': lambda x: format(x, '8.5E')})

    dat = dat*0.00101971621 # convert to G

    if fname[1] == 'EW2' or 'EW1': 
       s = 'EW' 
    elif fname[1] == 'NS2' or 'NS1':
       s = 'NS'
    else:
        s = 'HCOMB'
    #does filename already exist? If so delete and make new file.
    checkFileExist(fname[0]+s+'.AT2')
    #create the file for writing-fill with useless text 
    with open(fname[0]+s+'.AT2', 'w') as f:
        f.write('PEER NGA STRONG MOTION DATABASE RECORD\n')
        f.write('STUFF\n')
        f.write(
        'ACCELERATION TIME HISTORY IN UNITS OF G\n')
        f.write(str(npts)+'    '+str(dt)+'    NPTS, DT\n')    
        [f.write(str(s).strip('[ ]')+'\n') for s in dat]
    print('Conversion complete.')
  
def read_kiknet(fname):
    """Takes a kiknet file as input (e.g. foo.EW1) and returns the data 
    (in counts-float), the magnitude (jma-float), the station name (string), 
    the scaling factor (cm/s/s) and sampling rate (dt - [s]). 
    The data is returned as a N*8 matrix where N is the length of the matrix
    [e.g. len(data)]. The number of samples = N*8"""
   
    #open the file - check if gzipped or not - open appropriately
    if fname[-3:] = '.gz'
        f = gzip.open(fname, 'rt')
    else:
        f = open(fname, 'rt')
    strings = f.readlines() # reads each line of file into a list of the lines as text
    f.close() #close the file
    #removes '\n' and splits each line of list into separate strings
    strings = [s.strip() for s in strings]
    #station name
    stname = strings[5][2]
    #magnitude (jma)
    jmamag = float(strings[4][1])
    #frequency/dt
    freq = strings[10][2].strip('Hz')
    dt = 1/float(freq)
    #scaling factor block
    dump = strings[13][2].split('(gal)/') #ARRGH WHY FORMAT LIKE THIS
    scalingF = float(dump[0])/float(dump[1])
    #origin date and time
    odate_time = strings[0][2]+' '+strings[0][3]
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

    return dat, stname, jmamag, dt, scalingF, odate_time

#combine horz components and apply calibration 
def combine_comps(datEW, datNS):
    """Combines the two time series - magnitude of vectors. """

    rows = len(datEW)

    datEW = sg.detrend(datEW.reshape(1, rows*8)[0])
    datNS = sg.detrend(datNS.reshape(1, rows*8)[0])
 
    if np.min(datEW) < np.min(datNS):
        lowest = int(np.min(datEW))
    else:
        lowest = int(np.min(datNS))

    datEW = datEW - (lowest-10)
    datNS = datNS - (lowest-10)

    comb = (datEW**2 + datNS**2)**(1/2)
    
    return comb.reshape(rows, 8)-comb[0]

def check4BasicInput():
    """Checks for the minimum input of two argument variables for the script to
       operate normally."""
    try: 
        fname = str(sys.argv[1]) 
        fname = str(sys.argv[2]) 
   
    except IndexError:
        print(
        'Please specify input operation and file name. e.g. \'-p foo1234.EW1\'')
        print('Type: $python conv4strata.py --help for basic instructions')
        print('Program terminated.')
        sys.exit()

def check4AdvancedInput():
    """Checks for arguments required for more advanced operations to 
       function."""
    try: 
        fname = str(sys.argv[3]) 
         
   
    except IndexError:
        print(
        'Please specify second file name. e.g. \'-p foo1234.NS1\'')
        print('Program terminated.')
        sys.exit()    

def checkFileExist(fname):
    if os.path.isfile(fname):
        print('File already exists: File Deleted.')
        os.remove(fname)
    else:
        print('Filename does not exist: Resuming.')
    
def check4AnyInput():
    if len(sys.argv) <= 1:
        print('No inputs given.Type: $python conv4strata.py --help for basic instructions')
        sys.exit()  

def check4CorrectInput():
    options = ['--help', '-p', '-hp', '-c','-hc']
   
    if sum([s==sys.argv[1] for s in options]) == 0:
        print('Invalid Option.Type: $python conv4strata.py --help for basic instructions')
        sys.exit()      
#---------------------------------run cmds-------------------------------------#
#help module
#if str(sys.argv[1] == '--help'):
#print('Help module for conv4strata')

#preview waveform
if __name__ == "__main__":


    check4AnyInput()
    check4CorrectInput()  
 
    if str(sys.argv[1]) == '--help':
        print('$ python conv4strata -p name_of_file.EW2 to preview waveform')
        print(
        '$ python conv4strata -hc name_of_file.EW2 name_of_file.NS2 to combine waveforms and preview.')
        print(
        '$ python conv4strata -c name_of_file.EW2 to convert to .AT2. Outputs to current directory.')
        print(
        '$ python conv4strata -hc name_of_file.EW2 name_of_file.NS2 to combine waveforms and convert to .AT2. Outputs to current directory.')
        
    if str(sys.argv[1]) == '-p':
        check4BasicInput()
        preview_plot()
        
    if str(sys.argv[1]) == '-hp':
        check4AdvancedInput()
        preview_horz()
    
    if str(sys.argv[1]) == '-c':
        check4BasicInput()
        simple_conv() 
    
    if str(sys.argv[1]) == '-hc':
        check4AdvancedInput()             
        horz_conv()
    

#if str(sys.argv[1]) == '-h':
   #EWfile = sys.argv[2]
   #NSfile = sys.argv[3]
   #Note - obspy gives values in m/s/s after apply calib const
   #st = read(EWfile)
   #st += read(NSfile)
   #st.detrend()
   #combine the horizontal componants
   #horz = combine_comps(st[0].data, st[1].data)*st[0].stats.calib*100
   #detrend the combined signal
   #horz = sg.detrend(horz)


#if str(sys.argv[1]) == '-s':
   #f = sys.argv[2] 





        

