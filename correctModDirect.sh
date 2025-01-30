#!/bin/bash
# Correct Filament Model Direction
# Usage: correctModDirect.sh modelFile direction correctedModelFile

module load imod

#echo $1
#echo $2
#echo $3

basename=${1%.*}

model2point -Contour $1 ${basename}.txt

if [ $2 -eq 0 ]
then
  echo "sort -k1,1n -k 3n ${basename}.txt > ${basename}_fix.txt"
  sort -k1,1n -k 3n ${basename}.txt > ${basename}_fix.txt
else
  echo "sort -k1,1n -k 3nr ${basename}.txt > ${basename}_fix.txt"
  sort -k1,1n -k 3nr ${basename}.txt > ${basename}_fix.txt
fi

point2model ${basename}_fix.txt $3
