#!/bin/bash
# Script to count particles in particles folder

# Threshold for display
maxno=15

echo "Display folder with less than $maxno particles"

for folder in particles/*
do
	count=`ls $folder/particle* | wc -l`
	if [ $count -lt $maxno ]
	then
		echo $folder $count 
	fi
done
