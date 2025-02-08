
#!/bin/bash
# createMrcTxt.sh [path_to_models] [regular_expression_of_mrc_files] [project_Path]
# Example: createMrcTxt.sh tomograms "*_MT.mrc" .

echo "$1"
echo "$2"
echo "$3"
cd $1
rm $3mrcfiles.txt
ls $2 >> $3mrcfiles.txt
