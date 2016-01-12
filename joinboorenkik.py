#SCRIPT TO CONCATENATE THE BOORE DISTANCE METRICS WITH KIK DATABASE


#IMPORT MODULES
import numpy as np
from sys import argv

#CALL FILENAME AS ARGUMENT VARIABLE TO LOAD IT INTO SCRIPT
script, fname, boore1, boore2, name = argv

def linecounter(fname):
    with open(fname) as f:
       lc = sum(1 for _ in f)
       name = f.name
    return lc, name

def extractEWNSUDonetwo(fname, lc):
    """ This function separates the file into EW1, NS1, UD1, EW2, NS2, UD2. 
    It takes the filename and the linecount of the file as arguments.  """
    a = np.loadtxt(fname, usecols=(1,2,3,4,5,6))
    #calculate number of rows there should be for each channel
    rows = (lc / 6)
    EWone = a[0:rows, 0:6]
    NSone = a[(rows):(rows * 2), 0:6]
    UDone = a[(rows * 2):(rows * 3), 0:6]
    EWtwo = a[(rows * 3):(rows * 4), 0:6]
    NStwo = a[(rows * 4):(rows * 5), 0:6] 
    UDtwo = a[(rows * 5):(rows * 6), 0:6] 
    return EWone, NSone, UDone, EWtwo, NStwo, UDtwo

def epiDistSort(data):
    """ This function sorts data in terms of epicentral distance which is 
 column 8 in the original file. """
    episorted = data[data[:, 4].argsort()]
    return episorted



def concatKikBoore(kik, boore):
    b = np.loadtxt(boore)
    
    c = np.concatenate((kik, b), axis=1)
    return c

def builddatabase(EW1, NS1, UD1, EW2, NS2, UD2, name):
    c = np.concatenate((EW1, NS1, UD1, EW2, NS2, UD2), axis=0)
    np.savetxt(name, c, fmt='%10.5f')




EWone, NSone, UDone, EWtwo, NStwo, UDtwo = extractEWNSUDonetwo(fname, 
    (linecounter(fname)[0]))

EWoneEpi = epiDistSort(EWone)
NSoneEpi = epiDistSort(NSone)
UDoneEpi = epiDistSort(UDone)
EWtwoEpi = epiDistSort(EWtwo)
NStwoEpi = epiDistSort(NStwo)
UDtwoEpi = epiDistSort(UDtwo)

EW1 = concatKikBoore(EWoneEpi, boore1)
NS1 = concatKikBoore(NSoneEpi, boore1)
UD1 = concatKikBoore(UDoneEpi, boore1)

EW2 = concatKikBoore(EWtwoEpi, boore2)
NS2 = concatKikBoore(NStwoEpi, boore2)
UD2 = concatKikBoore(UDtwoEpi, boore2)

builddatabase(EW1, NS1, UD1, EW2, NS2, UD2, name)
