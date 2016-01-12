#SCRIPT TO SEPATATE/SORT THE INITIAL DATABASE AND SOME PLOTTING FEATURES

#IMPORT MODULES
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
import scipy.signal as sg


script, fname = argv

def linecounter(fname):
    with open(fname) as f:
       lc = sum(1 for _ in f)
       name = f.name
    return lc, name

def extractEWNSUDonetwo(f, lc):
    """ This function separates the file into EW1, NS1, UD1, EW2, NS2, UD2. 
    It takes the filename and the linecount of the file as arguments.  """
    a = np.loadtxt(f, usecols=(1,2,3,4,5,6))
    #calculate number of rows there should be for each channel
    rows = (lc / 6)
    EWone = a[0:rows, 0:6]
    NSone = a[(rows + 1):(rows * 2), 0:6]
    UDone = a[(rows * 2 + 1):(rows * 3), 0:6]
    EWtwo = a[(rows * 3 + 1):(rows * 4), 0:6]
    NStwo = a[(rows * 4 + 1):(rows * 5), 0:6] 
    UDtwo = a[(rows * 5 + 1):(rows * 6), 0:6] 
    return EWone, NSone, UDone, EWtwo, NStwo, UDtwo



def epiDistSort(data):
    """ This function sorts data in terms of epicentral distance which is 
 column 8 in the original file. """
    episorted = data[data[:, 4].argsort()]
    return episorted

def plotmaxAvepiD(EW, NS, UD, xMin, xMax):
    plt.subplot(3,1,1).set_xlim([xMin, xMax])
    plt.plot(EW[:,4], EW[:,3], '.') #6th column is max acc as a fnctn of epiD
    plt.subplot(3,1,2).set_xlim([xMin, xMax])
    plt.plot(NS[:,4], NS[:,3], '.')
    plt.subplot(3,1,3).set_xlim([xMin, xMax])
    plt.plot(UD[:,4], UD[:,3], '.')
    plt.show()

def plotmaxAvepiDloglog(EW, NS, UD, xMin, xMax):
    plt.subplot(3,1,1).set_xlim([xMin, xMax]) 
    plt.loglog(EW[:,4], EW[:,3], '.') #6th column is max acc as a fnctn of epiD
    plt.subplot(3,1,2).set_xlim([xMin, xMax]) #8th column is epicentral dist
    plt.loglog(NS[:,4], NS[:,3], '.')
    plt.subplot(3,1,3).set_xlim([xMin, xMax]) 
    plt.loglog(UD[:,4], UD[:,3], '.')
    plt.show()



   
###############################################################################
#TESTING SECTION
# Extract the different directions from the database 
EWone, NSone, UDone, EWtwo, NStwo, UDtwo = extractEWNSUDonetwo(fname, 
    (linecounter(fname)[0]))

# Sort by epicentral distance 
EWoneEpi = epiDistSort(EWone)
NSoneEpi = epiDistSort(NSone)
UDoneEpi = epiDistSort(UDone)
EWtwoEpi = epiDistSort(EWtwo)
NStwoEpi = epiDistSort(NStwo)
UDtwoEpi = epiDistSort(UDtwo)

np.savetxt('EWone2boore.txt', EWone[0:(linecounter(fname)[0]/6), 0:3], fmt='%10.5f')

np.savetxt('EWtwo2boore.txt', EWtwo[0:(linecounter(fname)[0]/6), 0:3], fmt='%10.5f')
#plotmaxAvepiD(EWoneEpi, NSoneEpi, UDoneEpi, 19, 110)

#plotmaxAvepiDloglog(EWoneEpi, NSoneEpi, UDoneEpi, 10, 110)




