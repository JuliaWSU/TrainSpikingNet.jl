#!/bin/bash

# Set the filename
filename='test.jld2'
# Create an empty file
touch $filename
# Check the file is exists or not
if [ -f $filename ]; then
   rm test.txt
   echo "$filename is removed"
fi
#if [ -f 'test.jld2']; then
#    rm ~/data_ben/*.jld2
#fi
./tsn.sh init ~/data_ben
./tsn.sh train gpu ~/data_ben
./tsn.sh test gpu ~/data_ben
./tsn.sh plot ~/data_ben/test.jld2
