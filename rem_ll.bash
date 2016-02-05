#!/bin/bash

#WARNING, ONLY RUN THIS SCRIPT ONCE TO AVOID LOSING DATA UNECCESSARILY

for file in *.*[0-9]; 
    do echo $file; 
    sed -i '' -e '$ d' $file; #for mac
done
