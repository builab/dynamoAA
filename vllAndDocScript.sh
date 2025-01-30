#!/bin/bash
# Updated by Mike Strauss 2022/06/08

args=("$@")

index=1

rm $5
rm $6

firstFound=false
first=1
last=2

cat $3 | while read line || [[ -n $line ]];
do
    a=$(dirname $line)
    #echo "$a"
    # cp $1$a/$a$7 $2cat $6 | while read line || [[ -n $line ]];
    tiltfile=$1$a/$a.tlt  # this defines the tilt file
    testfirst=$(head -1 $tiltfile)
    testlast=$(tail -1 $tiltfile)

    # test if the values were found
    if [ -z "$testfirst" ] || [ -z "$testlast" ] ;
        then
        echo "tilt range not found, cannot find tilts in $tiltfile.  Using range: $first, $last"
    else
        # set values only if they exist
        first=$testfirst
        last=$testlast
        echo $a $first $last
    fi
    
    #while read this;
    #do
    #    # echo "$this"
    #    if [[ "$this" == TiltAngle* ]];
    #    then
    #        if [[ $firstFound == false ]];
    #        then
    #            first=${this:12}
    #            # echo $first
    #            firstFound=true
    #        fi
    #        last=${this:12}
    #        # echo $last
    #    fi
    #done < $1$a/$a.mrc.mdoc
    


    VARNAME="${first%%[[:cntrl:]]}"
    echo "$index $1$a/$a$7" >> $5
    echo "$1$a/$a$7" >> $6
    echo "* ytilt =" $VARNAME $last >> $6
    echo "* apix =" $8 >> $6
    index=$((index + 1))
done
