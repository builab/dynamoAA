#!/bin/bash
# createModTxt.sh [path_to_models] [regular_expression_of_mod_files] [project_Path]
# createModTxt.sh tomograms "*_MT.mod" .

echo "$1"
echo "$2"
echo "$3"
cd $1
rm $3modfiles.txt
ls $2 >> $3modfiles.txt
