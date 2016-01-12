#!/bin/bash

#CHANGE OUTPUT.txt TO SUITABLE FILE NAME AND CHANGE THE DIP ANGLE AND MAGNITUDE

for file in *.EW1 *.NS1 *.UD1; do
    echo "$file $(python /Users/jamesholt/seismograms/process_kikhdr.py $file.h 90 6.8)" >> ./M7_2000_10_06_tottoribest.txt
done

for file in *.EW2 *.NS2 *.UD2; do
    echo "$file $(python /Users/jamesholt/seismograms/process_kikhdr.py $file.h 90 6.8)" >> ./M7_2000_10_06_tottoribest.txt
done
