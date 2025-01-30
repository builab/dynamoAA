#!/bin/bash
# model2pointscript.sh modDir modelDir modelfile stringToBeRemoved

module load imod

args=("$@")

cat $3 | while read line || [[ -n $line ]];
do
    a=$(dirname $line)
    echo "$a"
    suffix=$4
    echo "$suffix"
    b=$(echo "$line" | sed -e "s/$suffix$//")
    echo "$b"
    model2point -Contour $1$line $2$b.txt
    echo "model2point -Contour $1$line $2$b.txt"
done
