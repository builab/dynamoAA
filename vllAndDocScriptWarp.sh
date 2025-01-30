#!/bin/bash
# Updated by Mike Strauss 2022/06/08
# Important: Copy all the tlt file in the same directory as the rec & mod
# vllAndDocScriptWarp.sh [models_dir] [listOfTomograms_dir] [modelfile.txt] stringToBeRemoved docFile vllFile recSuffix apixel
# vllAndDocScriptWarp.sh tomograms/ catalogs/list_of_tomograms modfiles.txt _8.48Apx_MT.mod tomograms.doc tomograms.vll _8.48Apx.rec 8.48
# Warp has revsere angles

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
    b=$(basename $line)
    b=${b%$4}

    #echo "$a"
    # cp $1$a/$a$7 $2cat $6 | while read line || [[ -n $line ]];
    tiltfile=$1$a/$b.tlt  # this defines the tilt file
    testfirst=$(head -1 $tiltfile)
    testlast=$(tail -1 $tiltfile)

    # test if the values were found
    if [ -z "$testfirst" ] || [ -z "$testlast" ] ;
        then
        echo "tilt range not found, cannot find tilts in $tiltfile.  Using range: $first, $last"
    else
        # set values only if they exist
        # Reverse value due to Warp different
	first=$(echo "$testfirst * -1" | bc)
	last=$(echo "$testlast * -1" | bc)
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
    echo "$index $1$a/$b$7" >> $5
    echo "$1$a/$b$7" >> $6
    echo "* ytilt =" $VARNAME $last >> $6
    echo "* apix =" $8 >> $6
    index=$((index + 1))
done
