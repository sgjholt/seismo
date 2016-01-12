#!/bin/bash

for file in *.*[0-9]; do
    echo "$file $(python plot_waveforms.py $file $file.h)" 
done
