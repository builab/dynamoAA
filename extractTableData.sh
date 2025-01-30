#!/bin/bash

args=("$@")

if [ -f "$1" ] ; then
    rm "$1"
fi

if [ -f "$2" ] ; then
    rm "$2"
fi

x="$3$4_1/aligned.tbl"
echo "Loading $4"

cat $x | while read line || [[ -n $line ]];
do
    a=$(echo $line | awk '{print $4+$24}')
    b=$(echo $line | awk '{print $5+$25}')
    c=$(echo $line | awk '{print $6+$26}')
    #echo "$a, $b, $c"
    echo "$a, $b, $c" >> $1
    
done

y="$3$4_2/aligned.tbl"
cat $y | while read line || [[ -n $line ]];
do
    a=$(echo $line | awk '{print $4+$24}')
    b=$(echo $line | awk '{print $5+$25}')
    c=$(echo $line | awk '{print $6+$26}')
    #echo "$a, $b, $c"
    echo "$a, $b, $c" >> $2
    
done
