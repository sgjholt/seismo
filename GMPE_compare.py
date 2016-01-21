
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

# Tracking argvs - argv[1] = condition for coeff (R, SS or N)
#                - argv[2] =  condition for coeffs (crustal or slab)
#                - argv[3] = magnitude of earthquake (Mw)



def coeffPGApicker(condition1, condition2):
    a, b, c, d, e, Sr, Si, Ss, Ssl = [
    1.101, -0.00564, 0.0055, 1.080, 0.01412, 0.251, 0.000, 0.2607, -0.528]
    
    if condition1 == 'R' and condition2 == 'crustal': 
        Si, Ss, Ssl = [0, 0, 0]
        print('Selected coefficients for %s and %s type event.' %(
        condition1, condition2))
        return a, b, c, d, e, Sr, Si, Ss, Ssl
        
    elif condition1 == 'R' and condition2 == 'slab':
        print('Selected coefficients for %s and %s type event.' %(
        condition1, condition2))
        return a, b, c, d, e, Sr, Si, Ss, Ssl     
    
    elif condition1 == 'SS' or condition1 == 'N' and condition1 == 'crustal':
        Sr, Si, Ss, Ssl = [0, 0, 0, 0]
        print('Selected coefficients for %s and %s type event.' %(
        condition1, condition2))
        return a, b, c, d, e, Sr, Si, Ss, Ssl
    else:
        print(
        'Error. You have not selected valid conditions for coefficients.')


def zhao_06_GMPE():
    
    a, b, c, d, e, Sr, Si, Ss, Ssl = coeffPGApicker(str(argv[1]), str(argv[2]))
    
    splitdb = splitcompANDloc(np.loadtxt(
    '/Users/jamesholt/Documents/MRes/20001006133000/20001006133000.kik/M6.7_data_boore_complete.txt'))
    
    d, x = distSorter(splitdb[0], 'rrup')
 
    for i in x:
    r = x[i] + c * np.exp(d * argv[3])
    y[i] = (
    (a * argv[3]) + (b * x[i]) - np.log10(r) + Sr + Si + Ss + (Ssl 
    * np.log10(x[i])) + 0.6 + 0.723
    



def splitcompANDloc(data):
    """ This function separates the file into EW1, NS1, UD1, EW2, NS2, UD2. 
    It takes the file as an argument.  """
    
    #calculate number of rows there should be for each channel
    rows = len(data) / 6
    EWone = data[0:rows, 0:(data.size/len(data))]
    NSone = data[(rows):(rows * 2), 0:(data.size/len(data))]
    UDone = data[(rows * 2):(rows * 3), 0:(data.size/len(data))]
    EWtwo = data[(rows * 3):(rows * 4), 0:(data.size/len(data))]
    NStwo = data[(rows * 4):(rows * 5), 0:(data.size/len(data))] 
    UDtwo = data[(rows * 5):(rows * 6), 0:(data.size/len(data))] 
    
    listit = [EWone, NSone, UDone, EWtwo, NStwo, UDtwo]
    return listit


def distSorter(data, sortingcol):
    if sortingcol == 'epi':
        col = 4
    elif sortingcol == 'hyp':
        col = 5
    elif sortingcol == 'rrup':
        col = 6
    elif sortingcol == 'rcmbl':
        col = 7
    elif sortingcol == 'rjb':
        col = 8
    else:
        print('Incorrect sorting column, please choose rrup, rcmbl or rjb.')
    
    b = np.array([data[3], data[col]])
    new = b.transpose()
    distsorted = new[new[:, 1].argsort()]
    
    x = distmaker(distsorted[1])
    return distsorted, x

def distmaker(data):
    x = np.arange(0, data.max())
    return x


    
    
    
        
    
