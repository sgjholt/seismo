#!/bin/bash
for file in *.*[0-9]; 
do echo $file; 
	
        awk 'BEGIN{}\
        {if($1 == "Origin" && $2 == "Time"){Odate=$3;Otime=$4;}\
        if($1 == "Lat."){lat=$2;}\
	    if($1 == "Long."){lon=$2;}\
        if($1 == "Depth." && $2 == "(km)"){depth=$3;}\
        if($1 == "Mag."){mag=$2}\
        if($1 == "Station" && $2 == "Code"){stncd=$3;}\
	    if($1 == "Station" && $2 == "Lat." ){stnlat=$3;}\
	    if($1 == "Station" && $2 == "Long."){stnlon=$3;}\
        if($1 == "Station" && $2 == "Height(m)"){stnht=$3;}\
        if($1 == "Record" && $2 == "Time"){recD=$3;recT=$4;}\
        if($1 == "Sampling" && $2 == "Freq(Hz)"){spF=$3;}\
        if($1 == "Duration" && $2 == "Time(s)"){durT=$3;}\
        if($1 == "Scale" && $2 == "Factor"){scaleF=$3;}\
        if($1 == "Max." && $2 == "Acc." && $3 == "(gal)"){maxAc=$4;}\
        };END{print Odate, Otime, lat, lon, depth, mag, stnlat, stnlon, stnht, recD, recT, spF durT, scaleF, maxAc}' $file > ${file}.h; 

sed 's/(gal)\// /g' ${file}.h > tmp.dat; mv tmp.dat ${file}.h
sed 's/Hz/ /g' ${file}.h > tmp.dat; mv tmp.dat ${file}.h
sed 's/\// /g' ${file}.h > tmp.dat; mv tmp.dat ${file}.h
sed 's/:/ /g' ${file}.h > tmp.dat; mv tmp.dat ${file}.h
done


