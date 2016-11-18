# Vs_variable- A SCRIPT TO CALCULATE Vs30/50/etc for the KiK-Net sites.
# VERSION: 0.9
# AUTHOR(S): JAMES HOLT - UNIVERSITY OF LIVERPOOL
#            BEN EDWARDS - UNIVERSITY OF LIVERPOOL 
# EMAIL: j.holt@liverpool.ac.uk; Ben.Edwards@liverpool.ac.uk
# LAST EDIT: 2016/11/11
# 
# 2016/11/07: BAREBONES
# 2016/11/11: Completed parsing function to read .dat files
#             Completed Vs_Variable function
#             Added ability to print answer to command line
#---------------------------------modules--------------------------------------#
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import csv
#import gzip


#---------------------------------core fns-------------------------------------#
def Vs_variable(fname,refD):

    dat = readKikVelFile(fname)
    #Check to see if there is a depth of refD (special case)

    #if refD > data[:,2][-1]:
      #print('Cannot choose reference depth deeper than available data.')
      #sys.exit()

    if sum([n == refD for n in dat[:,2]]) == 0:
 
        VsX = ((refD - np.sum(dat[:,1][dat[:,2]<refD])) + np.sum(
        dat[:,1][dat[:,2]<refD])) / (np.sum(
        dat[:,1][dat[:,2]<refD]/dat[:,4][dat[:,2]<refD]) + ((30 - np.sum(
        dat[:,1][dat[:,2]<refD]))/dat[:,4][dat[:,2]>refD][0]))
   
   #In the special case the indexing is easier
    if sum([n == refD for n in dat[:,2]]) == 1: 
        VsX = np.sum(dat[:,1][dat[:,2]<=refD]) / (np.sum(
        dat[:,1][dat[:,2]<=refD]/dat[:,4][dat[:,2]<=refD])) 

    return np.array([refD,VsX])    

def AllSites(refD,option):
    path = '/data/share/Japan/SiteInfo/'
    fnames = sorted(glob.glob(path+'physicalData/*.dat'))

    if option == 'newfile':
        with open(path+'Vs'+str(refD)+'.csv', 'wt') as f:
            f.write('Site Code, Vs'+str(refD)+' (m/s)'+'\n')
    
        with open(path+'Vs'+str(refD)+'.csv', 'at') as f:
            [f.write(fname[0:6]+','+str(
            Vs_variable(fname, refD)[1])+'\n') for fname in fnames]

    if option == 'addtocsv':
        del refD

        with open(path+'SiteINFO.csv', 'r') as f:
            strings = [s.split(',') for s in [b.strip('\n') for b in f.readlines()]]


        with open(path+'SiteINFOext.csv', 'w') as f:
            writer = csv.writer(f, lineterminator='\n')
            for refD in [5,10,20,30,50,100,200]:
                strings[0].append('Vs'+str(refD))
                row = 1
                for fname in fnames:
                    strings[row].append(str(Vs_variable(fname, refD)[1]))
                    row += 1
            writer.writerows(strings)
#-----------------------------------utils--------------------------------------#

def readKikVelFile(fname):
    """Reads in a site velocity file in the SITE.dat format, returns 
    numpy array containing the data without the headers. Some files are 
    missing values in the first row and the bottom line always has '  -----,'
    for thickness depth. The function replaces blank lines with the velocity
    in the cells below and the thickness with 1000 m and depth with final 
    depth measurement + 1000 m. Returned is a square N*M matrix where M is 
    5 and N is arbitrary. Returned file has no headers - Col0 = Layer Num
    Col1 = Thickness (m), Col2 = Depth (m), Col3 = Vp (m/s) and Col4 = Vs (m/s) 
    """
    with open(fname) as f:
        strings = [s.split(',') for s in [b.strip('\n') for b in f.readlines()]]
    
    if strings[-1][1] == '   -----':
        strings[-1][1] = '1000'
        strings[-1][2] = str(float(strings[-2][2]) + 1000) 

    #sometimes true that you're missing values in top line (found 1 case so far)
    if strings[2][4] == '    0.00'  or strings[2][4] == '        ':
        strings[3][1] = strings [3][2] #make sure thickness is same as depth
        del strings[2]  #remove the whole line   
    #always true, the last thickness/depth will always look like this

    #make sure lines below also contain relevent data, if not remove the line
    indices = [i for i, s in enumerate(strings) if '        ' in s]
    if len(indices) == 1 and indices[0] == (len(strings)-1):
        del strings[indices[0]]



    #pre allocate numpy array and assign values to it     
    dat = np.zeros(((len(strings)-2),5)) 
    for n in range(2, len(strings)):
        dat[n-2,:]=[float(s) for s in strings[n]]

    return dat

def checkFiles():

    fnames = sorted(glob.glob('*.dat'))
    names = []
    for fname in fnames:
        with open(fname) as f:
            strings = [s.split(',') for s in [b.strip('\n') for b in f.readlines()]]
        if strings[2][4] ==  '    0.00' or strings[2][4] == '        ':
            strings[3][1] = strings [3][2] #make sure thickness is same as depth
            del strings[2]  #remove the whole line
        indices = [i for i, s in enumerate(strings) if '        ' in s]
        if len(indices) >= 1:
            print('Problem with {0}, consider editing.'.format(fname))


            names.append(fname)

    return names


#---------------------------------run cmds-------------------------------------#



if __name__ == "__main__":



    if str(sys.argv[1]) == '-p':
        fname = str(sys.argv[2])
        refD = int(sys.argv[3])
        print('Reference Depth (m), Vs'+str(refD)+ ' (m/s)')
        print(str(Vs_variable(fname,refD)[0])+',', str(Vs_variable(fname,refD)[1]))

    if str(sys.argv[1]) == '-all':
        print('Calculating Vs5,10,20,30,50,100,200 for all sites')
        print('Outputting to /SiteInfo/SiteINFOext.csv')
        AllSites('n/a', 'addtocsv')






