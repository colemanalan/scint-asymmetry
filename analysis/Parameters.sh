#!/bin/bash

#for ii in 15.0 15.5 16.0 16.5
#do 
	for jj in 0 10 20 30 40 50
	do
		#echo -n "$ii ;" ;
		echo -n "$jj ;"
		python3 Slope.py /home/ddazamarroquin/work/scint-asymmetry/data/new-i3-files/proton/lgE_16.0/Zen_$jj/000*.i3.gz
	done >Parameters_16_0.csv
#done > All_data.csv